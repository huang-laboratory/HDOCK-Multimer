import numpy as np
from scipy.spatial.transform import Rotation
from argparse import ArgumentParser
import os
import shutil
import copy
from datetime import datetime
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc
import tempfile
from contextlib import contextmanager
from pathlib import Path
import subprocess
import sys

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, "..", ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))

sys.path.append(BASE_PATH)
from parse_pdb import PDBParser, load_format_line, write_pdb_file


pdb_parser = PDBParser()
format_lines = load_format_line(os.path.join(BASE_PATH, "format.pdb"))
steric_clash_dist = 2
openmm_script = os.path.join(TOOL_PATH, "refine_model_openmm.py")
scorecom = os.path.join(TOOL_PATH, "scorecom.sh")
chdock = os.path.join(TOOL_PATH, "chdock")
compcn = os.path.join(TOOL_PATH, "compcn")
splitmodels_script = os.path.join(TOOL_PATH, "splitmodels.py")

current_dir = os.getcwd()
temp_dir = os.path.join(current_dir, "temp")
os.makedirs(temp_dir, exist_ok=True)


############    Basic Utility Functions    ############
# select appropriate number of parallel CPU cores
def get_num_cpus(max_workers: int):
    cpu_cap = int(os.getenv("NUM_CPUS", max_workers))
    cpu_count = os.cpu_count()
    return max(1, min(cpu_count, cpu_cap))

# renumber the chains in the complex starting from "A"
def change_chain_ids(initial_results):
    for result in initial_results:
        for i, chain in enumerate(result["chains"]):
            chain.id = chr(ord("A") + i)

# calculate rotation matrix for aligning two sets of coordinates using Kabsch's algorithm
def kabsch_rot(coords1, coords2):
    H = np.dot(coords1.T, coords2)
    U, S, VT = np.linalg.svd(H)
    rot = np.dot(VT.T, U.T)
    if np.linalg.det(rot) < 0.0:
        VT[2, :] *= -1
        rot = np.dot(VT.T, U.T)
    return rot

# calculate centroid for a chain
def calculate_centroid(chain):
    coords = [atom.coord for atom in chain.get_atoms()]
    return np.mean(coords, axis=0)

# shift center of the chain group to the center of mass
def center_chains(chains):
    coords = [atom.coord for chain in chains for atom in chain.get_atoms()]
    centered_centroid = np.mean(coords, axis=0)
    for chain in chains:
        for atom in chain.get_atoms():
            atom.coord -= centered_centroid
    return chains

# calculate the rotation and translation matrix for aligning two chains
def calculate_align_transformation(from_chain, to_chain):
    coords1 = np.array([atom.coord for atom in from_chain.get_atoms()])
    coords2 = np.array([atom.coord for atom in to_chain.get_atoms()])
    if coords1.shape != coords2.shape:
        raise ValueError("The sizes of aligned chains are different.")
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    centered_coords1 = coords1 - centroid1
    centered_coords2 = coords2 - centroid2
    rot = kabsch_rot(centered_coords2, centered_coords1)
    trans = centroid2 - np.dot(centroid1, rot)
    return rot, trans

# calculate the rotation matrix for aligning two axis directions
def calculate_rotation_matrix(from_axis, to_axis):
    if np.allclose(from_axis, to_axis, atol=1e-6):
        return np.eye(3)
    v = np.cross(from_axis, to_axis)
    s = np.linalg.norm(v)
    c = from_axis @ to_axis

    if np.isclose(s, 0, atol=1e-6) and np.isclose(c, -1, atol=1e-6):
        if not np.isclose(from_axis[0], 1.0, atol=1e-6):
            perpendicular_axis = np.array([1, 0, 0])
        else:
            perpendicular_axis = np.array([0, 1, 0])
        perp_vector = np.cross(from_axis, perpendicular_axis)
        perp_vector = perp_vector / np.linalg.norm(perp_vector)
        v_skew = np.array([
            [0, -perp_vector[2], perp_vector[1]],
            [perp_vector[2], 0, -perp_vector[0]],
            [-perp_vector[1], perp_vector[0], 0]
        ])
        return np.eye(3) + 2 * v_skew @ v_skew

    v_skew = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    rotation_matrix = np.eye(3) + v_skew + (v_skew @ v_skew) * ((1 - c) / (s ** 2))
    return rotation_matrix

# calculate CA rmsd between two chains
def calculate_ca_rmsd(chain1, chain2):
    coords1 = np.array(
        [atom.coord for atom in chain1.get_atoms() if atom.name == "CA"]
    )
    coords2 = np.array(
        [atom.coord for atom in chain2.get_atoms() if atom.name == "CA"]
    )
    diff = coords1 - coords2
    ca_rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return ca_rmsd

# calculate number of spatial conflicts in the chain group
def count_clashes(chains):
    chains_coords = [
        np.array([atom.coord for atom in chain.get_atoms()]) for chain in chains
    ]
    total_clash_count = 0
    for i in range(len(chains_coords)):
        for j in range(i + 1, len(chains_coords)):
            coords1, coords2 = chains_coords[i], chains_coords[j]
            distances = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=2)
            clash_count = np.sum(distances < steric_clash_dist)
            total_clash_count += clash_count
    return total_clash_count

# read C2 symmetric structure
def load_c2_symmetry_mode(pdb_file):
    model = pdb_parser.get_structure("", pdb_file)[0]
    chains = list(model.get_chains())
    if len(chains) != 2:
        raise ValueError("The number of chains in the C2 symmetry structure is incorrect.")
    return chains

# read C3 symmetric structure
def load_c3_symmetry_mode(pdb_file):
    model = pdb_parser.get_structure("", pdb_file)[0]
    chains = list(model.get_chains())
    if len(chains) != 3:
        raise ValueError("The number of chains in the C3 symmetry structure is incorrect.")
    chains = center_chains(chains)
    return chains

# extract geometric information from C3 symmetric structure
def analysis_c3_symmetry_mode(c3_chains):
    centroids = [calculate_centroid(chain) for chain in c3_chains]
    c3_axis = np.cross(centroids[2] - centroids[0], centroids[1] - centroids[0])
    c3_axis = c3_axis / np.linalg.norm(c3_axis)
    coords = [atom.coord for chain in c3_chains for atom in chain.get_atoms()]
    projections = [c3_axis @ coord for coord in coords]
    d = max(projections) - min(projections)
    sliding_range = np.ceil(2 * d / 3.0) * 3.0
    clash_count = count_clashes(c3_chains)
    return c3_axis, sliding_range, clash_count

# check whether the C2 symmetry is satisfied between two chains
def check_c2_symmetry_rmsd(chain1, chain2):
    coords1 = [atom.coord for atom in chain1.get_atoms()]
    coords2 = [atom.coord for atom in chain2.get_atoms()]
    if len(coords1) != len(coords2):
        raise ValueError("Lengths of C2 chains to be inspected are different.")
    coords = coords1 + coords2
    atm_idx = len(coords1) // 2
    c2_axis = np.cross(coords1[0] - coords2[0], coords1[atm_idx] - coords2[atm_idx])
    c2_axis /= np.linalg.norm(c2_axis)
    angle = np.pi
    rotation_matrix = Rotation.from_rotvec(angle * c2_axis).as_matrix()
    rotated_chain = chain2.copy()
    centroid = np.mean(coords, axis=0)
    for atom in rotated_chain.get_atoms():
        atom.coord = rotation_matrix @ (atom.coord - centroid) + centroid
    ca_rmsd = calculate_ca_rmsd(chain1, rotated_chain)
    return ca_rmsd

# check whether different chain pairs satisfy C2 symmetry
def check_tetrahedral_rmsd(new_chains):
    all_chains = []
    for chains in new_chains:
        centered_chains = [chain.copy() for chain in chains]
        centered_chains = center_chains(centered_chains)
        all_chains.append(centered_chains)
    ca_rmsds = [
        check_c2_symmetry_rmsd(all_chains[1][1], all_chains[3][1]),
        check_c2_symmetry_rmsd(all_chains[1][2], all_chains[2][2]),
        check_c2_symmetry_rmsd(all_chains[2][0], all_chains[3][0])
    ]
    return max(ca_rmsds)

@contextmanager
def managed_tempfile(suffix="", delete_ok: bool = True):
    with tempfile.NamedTemporaryFile(
            suffix=suffix, dir=temp_dir, delete=delete_ok
    ) as temp_file:
        try:
            yield temp_file.name
        finally:
            temp_file.close()

# calculate IT-score for the complex directly
def calculate_itscore(chains):
    with managed_tempfile(".pdb", delete_ok=False) as model_file:
        write_pdb_file(chains, model_file, format_lines)
        result = subprocess.run(f"bash {scorecom} {model_file}", shell=True, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        it_score = float(result.stdout.strip())
        return it_score, model_file

# refine structure using OpenMM energy minimization, then calculate IT-score
def calculate_itscore_md(chains, md_steps):
    with managed_tempfile(".pdb") as model_file, \
            managed_tempfile(".pdb", delete_ok=False) as md_file:
        write_pdb_file(chains, model_file, format_lines)
        subprocess.run(f"python {openmm_script} {model_file} {md_steps} {md_file}",
                       shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = subprocess.run(f"bash {scorecom} {md_file}", shell=True, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        it_score = float(result.stdout.strip())
        return it_score, md_file

# refine the final output model, not used
def refine_output_model(i, result, model_file, md_steps):
    if md_steps > 0:
        file1 = f"{current_dir}/0-Tsym_{i + 1}.pdb"
        write_pdb_file(result["chains"], file1, format_lines)
        subprocess.run(f"python {openmm_script} {file1} {md_steps} {model_file}",
                       shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        os.remove(file1)
    else:
        write_pdb_file(result["chains"], model_file, format_lines)

# print the final output models, not used
def print_results(results, output_file, num_cpus, md_steps):
    stem, suffix = output_file.rsplit(".", 1)
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(
                refine_output_model, i, result, f"{stem}_{i + 1}.{suffix}", md_steps
            ) for i, result in enumerate(results)
        ]
        for future in as_completed(futures):
            future.result()
    gc.collect()


############    Functions for Constructing T-symmetric Systems     ############
############    Build T-symmetry using C2-C3 Combinations
# check whether specific C2 and C3 symmetric mode can form T-symmetry
def generate_new_chains(c2_chains, c3_chains, c2_idx, c3_idx, rmsd_threshold=5.0):
    rot_trans = [calculate_align_transformation(c2_chains[0], chain) for chain in c3_chains]
    anchor_chains = []
    for rot, trans in rot_trans:
        anchor_chain = c2_chains[1].copy()
        for atom in anchor_chain.get_atoms():
            atom.transform(rot, trans)
        anchor_chains.append(anchor_chain)
    new_chains = [c3_chains.copy()]
    for k in range(3):
        rot, trans = calculate_align_transformation(c3_chains[k], anchor_chains[k])
        expanded_chains = []
        for chain in c3_chains:
            expanded_chain = chain.copy()
            for atom in expanded_chain.get_atoms():
                atom.transform(rot, trans)
            expanded_chains.append(expanded_chain)
        new_chains.append(expanded_chains)
    ca_rmsd = check_tetrahedral_rmsd(new_chains)
    tetrahedral = ca_rmsd < rmsd_threshold
    if not tetrahedral:
        return None

    c3_axes = []
    all_chains = []
    for chains in new_chains:
        centroids = [calculate_centroid(chain) for chain in chains]
        c3_axis = np.cross(centroids[2] - centroids[0], centroids[1] - centroids[0])
        c3_axis = c3_axis / np.linalg.norm(c3_axis)
        c3_axes.append(c3_axis)
        all_chains.extend(chains)
    result = {
        "c2_idx": c2_idx, "c3_idx": c3_idx, "c3_axes": c3_axes,
        "chains": all_chains, "ca_rmsd": ca_rmsd
    }
    return result

# traverse all C2-C3 combinations and retain those can form T-symmetry
def generate_initial_results_1(c2_modes_all, c3_modes_all, num_cpus):
    initial_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(
                generate_new_chains, c2_mode, c3_mode, i, j, 3.0
            )
            for i, c2_mode in enumerate(c2_modes_all)
            for j, c3_mode in enumerate(c3_modes_all)
        ]
        for future in as_completed(futures):
            result = future.result()
            if result:
                initial_results.append(result)
    gc.collect()

    if initial_results:
        change_chain_ids(initial_results)
    return initial_results

# perform sliding optimization on the C2-C3-formed T-symmetric systems
def sliding_optimization_1(result, md_steps):
    c3_axes = result["c3_axes"]
    start_chains = result["chains"]
    step_size = 1.0
    num_steps = int(5.0 / step_size)

    # stime = datetime.now()
    sliding_results = []
    for step in range(-num_steps, num_steps + 1):
        sliding_step = step * step_size
        new_chains = []
        for i, c3_axis in enumerate(c3_axes):
            for chain in start_chains[3 * i: 3 * (i + 1)]:
                new_chain = chain.copy()
                for atom in new_chain.get_atoms():
                    atom.coord += sliding_step * c3_axis
                new_chains.append(new_chain)
        sliding_results.append({"step": step, "chains": new_chains})
    if md_steps > 0:
        for sliding_result in sliding_results:
            it_score, file = calculate_itscore_md(sliding_result["chains"], md_steps)
            sliding_result["it_score"] = it_score
            sliding_result["file"] = file
    else:
        for sliding_result in sliding_results:
            it_score, file = calculate_itscore(sliding_result["chains"])
            sliding_result["it_score"] = it_score
            sliding_result["file"] = file
    sliding_results.sort(key=lambda x: x["it_score"])
    result["chains"] = sliding_results[0]["chains"]
    result["it_score"] = sliding_results[0]["it_score"]
    result["file"] = sliding_results[0]["file"]
    # etime = datetime.now()
    # i, j = optimized_result["c2_idx"], optimized_result["c3_idx"]
    # print(f"Sliding optimization in round ({i+1}, {j+1}) cost {etime - stime}")
    return result

# optimize all C2-C3-formed T-symmetric systems
def refine_results_1(initial_results, num_cpus, md_steps, num_models):
    refined_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(sliding_optimization_1, result, md_steps)
            for result in initial_results
        ]
        for future in as_completed(futures):
            result = future.result()
            refined_results.append(result)
            del result, future
    gc.collect()

    refined_results.sort(key=lambda x: x["it_score"])
    refined_results = refined_results[:min(len(refined_results), num_models)]
    return refined_results


############    Build T-symmetry using C3 trimers and Standard T-axes
# use the standard T-axes as the directions of the C3 axes in the T-symmetric system
def generate_tetrahedral_c3_axes(c3_axis):
    v1 = c3_axis / np.linalg.norm(c3_axis)
    standard_axes = np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]
    ])
    ref_axis = standard_axes[0] / np.linalg.norm(standard_axes[0])
    rotation_matrix = calculate_rotation_matrix(ref_axis, v1)
    tetrahedral_axes = [v1]
    for i in range(1, 4):
        ref_axis = standard_axes[i]
        rotated_axis = rotation_matrix @ ref_axis
        rotated_axis = rotated_axis / np.linalg.norm(rotated_axis)
        tetrahedral_axes.append(rotated_axis)
    return tetrahedral_axes

# check whether the two chains can satisfy C2 symmetry after rotation
def check_rotated_angle(chain1, chain2, axis, angle, rmsd_threshold=5.0):
    rotated_chain = chain2.copy()
    rotation_matrix = Rotation.from_rotvec(angle * axis).as_matrix()
    for atom in rotated_chain.get_atoms():
        atom.coord = rotation_matrix @ atom.coord
    ca_rmsd = check_c2_symmetry_rmsd(chain1, rotated_chain)
    if ca_rmsd < rmsd_threshold:
        return angle
    return None

# check whether the specific angle sets can ensure all chain pairs satisfy C2 symmetry
def validate_angles(fixed_chains, moving_chains, axes, angles, rmsd_threshold=5.0):
    all_chains = fixed_chains.copy()
    for chains, axis, angle in zip(moving_chains, axes, angles):
        rotation_matrix = Rotation.from_rotvec(angle * axis).as_matrix()
        for chain in chains:
            rotated_chain = chain.copy()
            for atom in rotated_chain.get_atoms():
                atom.coord = rotation_matrix @ atom.coord
            all_chains.append(rotated_chain)

    c2_pairs = [(0, 3), (1, 7), (2, 11), (4, 10), (5, 8), (6, 9)]
    total_rmsd = 0.0
    for i, j in c2_pairs:
        ca_rmsd = check_c2_symmetry_rmsd(all_chains[i], all_chains[j])
        if ca_rmsd >= rmsd_threshold:
            return None
        total_rmsd += ca_rmsd
    average_rmsd = total_rmsd / 6.0
    valid_angle_result = {"angles": angles, "rmsd": average_rmsd, "chains": all_chains}
    return valid_angle_result

# use specific C3 trimer and standard T-axes to build T-symmetric system
def copy_c3_chains(c3_idx, c3_chains, angle_step):
    c3_axis, sliding_range, clash_count = analysis_c3_symmetry_mode(c3_chains)
    c3_axes = generate_tetrahedral_c3_axes(c3_axis)
    new_chains = []
    for i in range(4):
        c3_axis = c3_axes[i]
        rotation_matrix = calculate_rotation_matrix(c3_axes[0], c3_axis)
        for chain in c3_chains:
            new_chain = chain.copy()
            for atom in new_chain.get_atoms():
                atom.coord = rotation_matrix @ atom.coord
            new_chains.append(new_chain)

    fixed_chains = new_chains[0: 3]
    moving_chains = [new_chains[3*i: 3*(i+1)] for i in range(1, 4)]
    num_steps = int(2 * np.pi / angle_step)
    sampled_angles = [[] for _ in range(3)]
    for i in range(3):
        for step in range(num_steps):
            angle = check_rotated_angle(
                fixed_chains[i], moving_chains[i][i], c3_axes[i + 1], step * angle_step
            )
            if angle:
                sampled_angles[i].append(angle)
        if len(sampled_angles[i]) == 0:
            return None
    valid_combinations = []
    for angle1 in sampled_angles[0]:
        for angle2 in sampled_angles[1]:
            for angle3 in sampled_angles[2]:
                angle_result = validate_angles(
                    fixed_chains, moving_chains, c3_axes[1:], [angle1, angle2, angle3]
                )
                if angle_result:
                    valid_combinations.append(angle_result)

    if valid_combinations:
        valid_combinations.sort(key=lambda x: x["rmsd"])
        all_chains = valid_combinations[0]["chains"]
        initial_result = {
            "c3_idx": c3_idx, "chains": all_chains, "c3_axes": c3_axes,
            "sliding_range": sliding_range, "clash_count": clash_count
        }
        return initial_result
    else:
        return None

# traverse all C3 symmetric modes and retain those can form T-symmetry
def generate_initial_results_2(c3_modes_all, num_cpus, angle_step):
    initial_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(copy_c3_chains, i, c3_mode, angle_step)
            for i, c3_mode in enumerate(c3_modes_all)
        ]
        for future in as_completed(futures):
            result = future.result()
            if result:
                initial_results.append(result)
    gc.collect()

    if initial_results:
        change_chain_ids(initial_results)
    return initial_results

# perform sliding optimization on the C3-formed T-symmetric systems
def sliding_optimization_2(idx, result, md_steps):
    c3_idx = result["c3_idx"]
    c3_axes = result["c3_axes"]
    start_chains = result["chains"]
    sliding_range = result["sliding_range"]
    clash_count = result["clash_count"]
    step_size = 3.0
    num_steps = int(sliding_range / step_size)

    # stime = datetime.now()
    sliding_results = []
    for step in range(-num_steps, num_steps + 1):
        sliding_step = step * step_size
        new_chains = []
        for i, c3_axis in enumerate(c3_axes):
            for chain in start_chains[3 * i: 3 * (i + 1)]:
                new_chain = chain.copy()
                for atom in new_chain.get_atoms():
                    atom.coord += sliding_step * c3_axis
                new_chains.append(new_chain)
        total_clash_count = count_clashes(new_chains)
        sliding_results.append(
            {"c2_idx": 0, "c3_idx": c3_idx, "step": step,
             "chains": new_chains, "clash": total_clash_count}
        )
    contact_indices = []
    for i, sliding_result in enumerate(sliding_results):
        if sliding_result["clash"] > 4 * clash_count:
            contact_indices.append(i)
    sliding_results = sliding_results[contact_indices[0] - 1: contact_indices[-1] + 2]
    for sliding_result in sliding_results:
        it_score, file = calculate_itscore(sliding_result["chains"])
        sliding_result["it_score"] = it_score
        sliding_result["file"] = file
    sliding_results.sort(key=lambda x: x["it_score"])
    if md_steps > 0:
        it_score, file = calculate_itscore_md(sliding_results[0]["chains"], md_steps)
        sliding_results[0]["it_score"] = it_score
        sliding_results[0]["file"] = file
    # etime = datetime.now()
    # print(f"Sliding optimization in round {idx + 1} cost {etime - stime}")
    return sliding_results[0]

# optimize all C3-formed T-symmetric systems
def refine_results_2(initial_results, num_cpus, md_steps, num_models):
    refined_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(sliding_optimization_2, i, result, md_steps)
            for i, result in enumerate(initial_results)
        ]
        for future in as_completed(futures):
            result = future.result()
            refined_results.append(result)
            del result, future
    gc.collect()

    refined_results.sort(key=lambda x: x["it_score"])
    refined_results = refined_results[:min(len(refined_results), num_models)]
    return refined_results


############    MAIN     ############
def main():
    parser = ArgumentParser()
    parser.add_argument("-n2", "--n2", type=int, default=10)
    parser.add_argument("-n3", "--n3", type=int, default=20)
    parser.add_argument("-w", "--workers", type=int, default=20)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("-m", "--md", type=int, default=500)
    parser.add_argument("-n", "--nmax", type=int, default=10)

    args = parser.parse_args()
    num_c2 = args.n2
    num_c3 = args.n3
    max_workers = args.workers
    num_cpus = get_num_cpus(max_workers)
    output_file = args.output
    stem, suffix = output_file.rsplit(".", 1)
    md_steps = args.md
    num_models = args.nmax

    start_time = time.time()
    # print(f"[INFO] Start T symmetric docking", flush=True)

    # Local symmetric modes sampling
    parent_dir = os.path.dirname(current_dir)
    mono_file = "A.pdb"
    cn_out_files = ["A_c2.out", "A_c3.out"]
    cn_pdb_files = ["A2.pdb", "A3.pdb"]
    cn_log_files = ["A2.log", "A3.log"]
    cn_types = [2, 3]
    cn_nums = [num_c2, num_c3]
    for i, cn_out_file in enumerate(cn_out_files):
        cn_pdb_file, cn_log_file = cn_pdb_files[i], cn_log_files[i]
        cn_type, cn_num = cn_types[i], cn_nums[i]
        if os.path.exists(f"{parent_dir}/{cn_out_file}"):
            shutil.copy(f"{parent_dir}/{cn_out_file}", current_dir)
            subprocess.run(
                f"{compcn} {cn_out_file} {cn_pdb_file} -nmax {cn_num} -complex > /dev/null &&"
                f"python {splitmodels_script} {cn_pdb_file}", shell=True, check=True
            )
        else:
            subprocess.run(
                f"{chdock} {mono_file} {mono_file} -symmetry {cn_type} -out {cn_out_file} > {cn_log_file} &&"
                f"{compcn} {cn_out_file} {cn_pdb_file} -nmax {cn_num} -complex > /dev/null &&"
                f"python {splitmodels_script} {cn_pdb_file}", shell=True, check=True
            )
            shutil.copy(cn_out_file, parent_dir)

    # Tetrahedral symmetry construction
    c2_modes_all = []
    c3_modes_all = []
    for i in range(1, num_c2 + 1):
        pdb_file = f"A2_{i}.pdb"
        c2_mode = load_c2_symmetry_mode(pdb_file)
        c2_modes_all.append(c2_mode)
    for i in range(1, num_c3 + 1):
        pdb_file = f"A3_{i}.pdb"
        c3_mode = load_c3_symmetry_mode(pdb_file)
        c3_modes_all.append(c3_mode)

    initial_results = generate_initial_results_1(c2_modes_all, c3_modes_all, num_cpus)
    if initial_results:
        selected_results = []
        c3_results = {}
        for result in initial_results:
            c3_idx = result["c3_idx"]
            if c3_idx not in c3_results:
                c3_results[c3_idx] = []
            c3_results[c3_idx].append(result)
        for c3_idx, results in c3_results.items():
            if not results:
                continue
            results.sort(key=lambda x: x["c2_idx"])
            selected_results.append(results[0])
        refined_results = refine_results_1(selected_results, num_cpus, md_steps, num_models)
        end_time = time.time()
        with open(f"{current_dir}/aa2.log", 'w+') as fout:
            for i, result in enumerate(refined_results):
                j, k = result["c2_idx"], result["c3_idx"]
                it_score = result["it_score"]
                print(f"Model {i + 1}, ({j + 1}, {k + 1}), IT-score: {it_score}", file=fout)
                model_file = f"{current_dir}/{stem}_{i + 1}.{suffix}"
                shutil.copy(result["file"], model_file)
            print(f"\nTetrahedral symmetric docking cost {end_time - start_time}",
                  file=fout, flush=True)

    else:
        print("No valid C2-C3 combination, try to construct Tsymmetry manually.")
        angle_step = np.pi / 12
        new_initial_results = generate_initial_results_2(c3_modes_all, num_cpus, angle_step)
        if new_initial_results:
            refined_results = refine_results_2(new_initial_results, num_cpus, md_steps, num_models)
            end_time = time.time()
            with open(f"{current_dir}/aa2.log", 'w+') as fout:
                for i, result in enumerate(refined_results):
                    j, k = result["c2_idx"], result["c3_idx"]
                    it_score = result["it_score"]
                    print(f"Model {i + 1}, ({j + 1}, {k + 1}), IT-score: {it_score}", file=fout)
                    model_file = f"{current_dir}/{stem}_{i + 1}.{suffix}"
                    shutil.copy(result["file"], model_file)
                print(f"\nTetrahedral symmetric docking cost {end_time - start_time}",
                      file=fout, flush=True)

        else:
            print("[WARN] No valid Tsymmetry configurations found.")
            print("[WARN] Try more C3 trimers or smaller angle steps.")
            end_time = time.time()

    shutil.rmtree(temp_dir)
    print(f"\nTetrahedral symmetric docking cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    main()

