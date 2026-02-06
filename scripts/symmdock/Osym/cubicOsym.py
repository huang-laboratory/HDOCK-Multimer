import sys
import numpy as np
from scipy.spatial.transform import Rotation
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
from argparse import ArgumentParser

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

# renumber the chains in the file starting from "A"
def change_file_chain_ids(input_file, output_file):
    chain_ids = ([chr(i) for i in range(ord("A"), ord("Z") + 1)] +
                 [chr(i) for i in range(ord("a"), ord("z") + 1)] +
                 [chr(i) for i in range(ord("0"), ord("9") + 1)])
    current_chain = None
    output_lines = []
    chain_idx = 0

    with open(input_file, 'r') as fin:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                if current_chain is None:
                    current_chain = chain_ids[chain_idx]
                new_line = line[:21] + current_chain + line[22:]
                output_lines.append(new_line)
            elif line.startswith("TER"):
                if current_chain is None:
                    current_chain = chain_ids[chain_idx]
                new_line = line[:21] + current_chain + line[22:]
                output_lines.append(new_line)

                chain_idx += 1
                if chain_idx > len(chain_ids):
                    print(f"[ERROR] Too many chains in input PDB or too few chain IDS provided.",
                          file=sys.stderr)
                    sys.exit(1)
                current_chain = chain_ids[chain_idx] if chain_idx < len(chain_ids) else None

    with open(output_file, 'w+') as fout:
        fout.writelines(output_lines)

# define the cyclic symmetric sampling process
def run_csym_dock(stem, idx, csym, cn_num):
    stem_file = f"{stem}_{idx}.pdb"
    out_file = f"{stem}_{idx}_c{csym}.out"
    log_file = f"cn_{idx}.log"
    out_pdb_file = f"{stem}_{idx}_{csym}.pdb"

    if not Path(stem_file).exists():
        print(f"[WARN] {stem_file} not found, skipping.", file=sys.stderr)
        return 0
    subprocess.run(
        f"{chdock} {stem_file} {stem_file} -symmetry {csym} -out {out_file} > {log_file} &&"
        f"{compcn} {out_file} {out_pdb_file} -nmax {cn_num} -complex > /dev/null &&"
        f"python {splitmodels_script} {out_pdb_file}", shell=True, check=True
    )
    os.remove(out_pdb_file)
    for j in range(1, 11):
        file1 = f"{stem}_{idx}_{csym}_{j}.pdb"
        if not Path(file1).is_file():
            raise RuntimeError(f"[ERROR] Change chain IDs failed for {file1}.")
        file2 = f"{stem}{csym}_{idx}_{j}.pdb"
        change_file_chain_ids(file1, file2)
        os.remove(file1)

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

# calculate the rotation and translation matrix for aligning two chain groups
def calculate_align_transformation(from_chains, to_chains):
    if isinstance(from_chains, list):
        coords1 = np.array([atom.coord for chain in from_chains for atom in chain.get_atoms()])
    else:
        coords1 = np.array([atom.coord for atom in from_chains.get_atoms()])
    if isinstance(to_chains, list):
        coords2 = np.array([atom.coord for chain in to_chains for atom in chain.get_atoms()])
    else:
        coords2 = np.array([atom.coord for atom in to_chains.get_atoms()])
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
def load_c4_symmetry_mode(pdb_file):
    model = pdb_parser.get_structure("", pdb_file)[0]
    chains = list(model.get_chains())
    if len(chains) != 4:
        raise ValueError("The number of chains in the C4 symmetry structure is incorrect.")
    chains = center_chains(chains)
    return chains

# read C2-C3 symmetric structure
def load_c23_symmetry_mode(pdb_file):
    model = pdb_parser.get_structure("", pdb_file)[0]
    chains = list(model.get_chains())
    if len(chains) != 6:
        raise ValueError("The number of chains in the C2-C3 symmetry structure is incorrect.")
    chains = center_chains(chains)
    return chains

# extract geometric information from C4 symmetric structure
def analysis_c4_symmetry_mode(c4_chains):
    centroids = [calculate_centroid(chain) for chain in c4_chains]
    c4_axis = np.cross(centroids[-1] - centroids[0], centroids[1] - centroids[0])
    c4_axis = c4_axis / np.linalg.norm(c4_axis)
    coords = [atom.coord for chain in c4_chains for atom in chain.get_atoms()]
    projections = [c4_axis @ coord for coord in coords]
    d1 = max(projections) - min(projections)
    distances = np.linalg.norm(coords - np.dot(coords, c4_axis)[:, None] * c4_axis, axis=1)
    d2 = 2 * max(distances)  # horizontal width, related to geometric properties of Osym
    sliding_range = np.ceil((2 * d1 + d2) / 3.0) * 3.0
    clash_count = count_clashes(c4_chains)
    return c4_axis, sliding_range, clash_count

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

# generate chain indices for checking C2 symmetry
def match_c4_pairs(c23_chains, c4_chains):
    c4_pairs_all = [
        [(0, 3), (2, 5), (4, 1)],
        [(3, 0), (5, 2), (1, 4)],
        [(0, 5), (2, 1), (4, 3)],
        [(5, 0), (1, 2), (3, 4)]
    ]  # anchor indices
    rot, trans = calculate_align_transformation(c23_chains[0], c4_chains[0])
    rotated_chain1 = c23_chains[3].copy()
    rotated_chain2 = c23_chains[5].copy()
    for atom in rotated_chain1.get_atoms():
        atom.transform(rot, trans)
    for atom in rotated_chain2.get_atoms():
        atom.transform(rot, trans)
    ca_rmsds = [
        calculate_ca_rmsd(rotated_chain1, c4_chains[1]),
        calculate_ca_rmsd(rotated_chain1, c4_chains[-1]),
        calculate_ca_rmsd(rotated_chain2, c4_chains[1]),
        calculate_ca_rmsd(rotated_chain2, c4_chains[-1])
    ]
    pair_idx = np.argmin(ca_rmsds)
    c4_pairs1 = c4_pairs_all[pair_idx]
    if pair_idx < 2:
        c4_pairs2 = [
            [(1, 1), (1, 4), (2, 1), (2, 4)],
            [(2, 0), (2, 3), (3, 0), (3, 3)],
            [(3, 2), (3, 5), (1, 2), (1, 5)]
        ]  # indices used for checking Tsymmetry
    else:
        c4_pairs2 = [
            [(1, 1), (1, 2), (3, 1), (3, 2)],
            [(1, 3), (1, 4), (2, 3), (2, 4)],
            [(3, 0), (3, 5), (2, 0), (2, 5)]
        ]  # indices used for checking Tsymmetry
    return c4_pairs1, c4_pairs2

# check whether different chain groups satisfy C4 symmetry
def check_c4_symmetry_rmsds(chains):
    centered_chains = [chain.copy() for chain in chains]
    centered_chains = center_chains(centered_chains)
    centroids = [calculate_centroid(chain) for chain in centered_chains]
    c4_axis = np.cross(centroids[-1] - centroids[0], centroids[1] - centroids[0])
    c4_axis = c4_axis / np.linalg.norm(c4_axis)
    angle = np.pi / 2
    ca_rmsds = []
    for i, chain in enumerate(centered_chains):
        chain1 = centered_chains[(i - 1) % 4]
        chain2 = centered_chains[(i + 1) % 4]
        rotated_chain = chain.copy()
        rotation_matrix = Rotation.from_rotvec(angle * c4_axis).as_matrix()
        for atom in rotated_chain.get_atoms():
            atom.coord = rotation_matrix @ atom.coord
        ca_rmsd = min(
            calculate_ca_rmsd(chain1, rotated_chain), calculate_ca_rmsd(chain2, rotated_chain)
        )
        ca_rmsds.append(ca_rmsd)
    return ca_rmsds

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
        file1 = f"0-Osym_{i + 1}.pdb"
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


############    Functions for Constructing O-symmetric Systems     ############
############    Build T-symmetry using C2-C3-C4 Combinations
# check whether specific C2-C3 and C4 symmetric mode can form O-symmetry
def generate_new_chains(c23_chains, c4_chains, c2_idx, c3_idx, c4_idx, rmsd_threshold=5.0):
    c4_pairs1, c4_pairs2 = match_c4_pairs(c23_chains, c4_chains)
    rot_trans = [
        calculate_align_transformation(
            [c4_chains[0], c4_chains[1]], [c23_chains[i], c23_chains[j]]
        ) for i, j in c4_pairs1
    ]
    anchor_chains = []
    for rot, trans in rot_trans:
        rotated_chains = [c4_chains[2].copy(), c4_chains[3].copy()]
        for rotated_chain in rotated_chains:
            for atom in rotated_chain.get_atoms():
                atom.transform(rot, trans)
        anchor_chains.append(rotated_chains)
    new_chains = [c23_chains.copy()]
    # new_chains = [copy.deepcopy(c23_chains)]
    for k in range(3):
        i, j = c4_pairs1[k]
        rot, trans = calculate_align_transformation(
            [c23_chains[i], c23_chains[j]], anchor_chains[k]
        )
        expanded_chains = []
        for chain in c23_chains:
            expanded_chain = chain.copy()
            for atom in expanded_chain.get_atoms():
                atom.transform(rot, trans)
            expanded_chains.append(expanded_chain)
        new_chains.append(expanded_chains)
    c4_chains1 = [new_chains[i][j] for i, j in c4_pairs2[0]]
    c4_chains2 = [new_chains[i][j] for i, j in c4_pairs2[1]]
    c4_chains3 = [new_chains[i][j] for i, j in c4_pairs2[2]]
    ca_rmsds1 = check_c4_symmetry_rmsds(c4_chains1)
    ca_rmsds2 = check_c4_symmetry_rmsds(c4_chains2)
    ca_rmsds3 = check_c4_symmetry_rmsds(c4_chains3)
    ca_rmsd = min(ca_rmsds1 + ca_rmsds2 + ca_rmsds3)
    octahedral = ca_rmsd < rmsd_threshold
    if not octahedral:
        return None

    c23_axes = []
    all_chains = []
    for chains in new_chains:
        centroids = [calculate_centroid(chain) for chain in chains]
        c23_axis = np.cross(centroids[2] - centroids[0], centroids[4] - centroids[0])
        c23_axis = c23_axis / np.linalg.norm(c23_axis)
        c23_axes.append(c23_axis)
        all_chains.extend(chains)
    result = {
        "c2_idx": c2_idx, "c3_idx": c3_idx, "c4_idx": c4_idx,
        "c23_axes": c23_axes, "chains": all_chains, "ca_rmsd": ca_rmsd
    }
    return result

# traverse all C2-C3-C4 combinations and retain those can form O-symmetry
def generate_initial_results_1(c23_modes_all, c4_modes_all, num_cpus):
    initial_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(
                generate_new_chains, c23_mode["chains"], c4_mode,
                c23_mode["c2_idx"], c23_mode["c3_idx"], i+1, 10.0
            )
            for c23_mode in c23_modes_all
            for i, c4_mode in enumerate(c4_modes_all)
        ]
        for future in as_completed(futures):
            result = future.result()
            if result:
                initial_results.append(result)
    gc.collect()

    if initial_results:
        change_chain_ids(initial_results)
    return initial_results

# perform sliding optimization on the C2-C3-C4-formed O-symmetric systems
def sliding_optimization_1(result, md_steps):
    c23_axes = result["c23_axes"]
    start_chains = result["chains"]
    step_size = 1.0
    num_steps = int(5.0 / step_size)

    # stime = datetime.now()
    sliding_results = []
    for step in range(-num_steps, num_steps + 1):
        sliding_step = step * step_size
        new_chains = []
        for i, c23_axis in enumerate(c23_axes):
            for chain in start_chains[6 * i: 6 * (i + 1)]:
                new_chain = chain.copy()
                for atom in new_chain.get_atoms():
                    atom.coord += sliding_step * c23_axis
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
    # i, j, k = result["c2_idx"], result["c3_idx"], result["c4_idx"]
    # print(f"Sliding optimization in round ({i}, {j})-{k} cost {etime - stime}")
    return result

# optimize all C2-C3-C4-formed O-symmetric systems
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


############    Build O-symmetry using C4 trimers and Standard O-axes
# use the standard O-axes as the directions of the C4 axes in the O-symmetric system
def generate_octahedral_c4_axes(c4_axis):
    v1 = c4_axis / np.linalg.norm(c4_axis)
    standard_axes = np.array([
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1]
    ])
    ref_axis = standard_axes[0] / np.linalg.norm(standard_axes[0])
    rotation_matrix = calculate_rotation_matrix(ref_axis, v1)
    octahedral_axes = [v1]
    for i in range(1, 6):
        ref_axis = standard_axes[i]
        rotated_axis = rotation_matrix @ ref_axis
        rotated_axis = rotated_axis / np.linalg.norm(rotated_axis)
        octahedral_axes.append(rotated_axis)
    return octahedral_axes

# check whether the two chains can satisfy C2 symmetry after rotation
def check_rotated_angle_1(chain1, chain2, axis, angle, rmsd_threshold=5.0):
    rotated_chain = chain2.copy()
    rotation_matrix = Rotation.from_rotvec(angle * axis).as_matrix()
    for atom in rotated_chain.get_atoms():
        atom.coord = rotation_matrix @ atom.coord
    ca_rmsd = check_c2_symmetry_rmsd(chain1, rotated_chain)
    if ca_rmsd < rmsd_threshold:
        return angle
    return None

# check whether the specific chain pairs in the two chain groups can
# satisfy C2 symmetry after rotation
def check_rotated_angle_2(fixed_chains, moving_chains, axis, angle, rmsd_threshold=5.0):
    rotation_matrix = Rotation.from_rotvec(angle * axis).as_matrix()
    rotated_chains = []
    ca_rmsds = []
    for i, chain in enumerate(moving_chains):
        rotated_chain = chain.copy()
        for atom in rotated_chain.get_atoms():
            atom.coord = rotation_matrix @ atom.coord
        ca_rmsd = check_c2_symmetry_rmsd(fixed_chains[i], rotated_chain)
        if ca_rmsd >= rmsd_threshold:
            return None
        rotated_chains.append(rotated_chain)
        ca_rmsds.append(ca_rmsd)
    ca_rmsd = max(ca_rmsds)
    angle_result = {"angle": angle, "rmsd": ca_rmsd, "chains": rotated_chains}
    return angle_result

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

    c2_pairs = [(0, 4), (1, 9), (2, 14), (3, 19), (5, 18), (7, 10), (8, 15), (13, 16)]
    total_rmsd = 0.0
    for i, j in c2_pairs:
        ca_rmsd = check_c2_symmetry_rmsd(all_chains[i], all_chains[j])
        if ca_rmsd >= rmsd_threshold:
            return None
        total_rmsd += ca_rmsd
    average_rmsd = total_rmsd / 8.0
    valid_angle_result = {"angles": angles, "rmsd": average_rmsd, "chains": all_chains}
    return valid_angle_result

# use specific C4 tetramer and standard O-axes to build O-symmetric system
def copy_c4_chains(c4_idx, c4_chains, angle_step):
    c4_axis, sliding_range, clash_count = analysis_c4_symmetry_mode(c4_chains)
    c4_axes = generate_octahedral_c4_axes(c4_axis)
    new_chains = []
    for i in range(6):
        c4_axis = c4_axes[i]
        rotation_matrix = calculate_rotation_matrix(c4_axes[0], c4_axis)
        for chain in c4_chains:
            new_chain = chain.copy()
            for atom in new_chain.get_atoms():
                atom.coord = rotation_matrix @ atom.coord
            new_chains.append(new_chain)

    fixed_chains = new_chains[0: 4]
    moving_chains = [new_chains[4*i: 4*(i+1)] for i in range(1, 5)]
    num_steps = int(2 * np.pi / angle_step)
    sampled_angles = [[] for _ in range(4)]
    for i in range(4):
        for step in range(num_steps):
            angle = check_rotated_angle_1(
                fixed_chains[i], moving_chains[i][i], c4_axes[i+1], step * angle_step
            )
            if angle:
                sampled_angles[i].append(angle)
        if len(sampled_angles[i]) == 0:
            return None
    valid_combinations = []
    for angle1 in sampled_angles[0]:
        for angle2 in sampled_angles[1]:
            for angle3 in sampled_angles[2]:
                for angle4 in sampled_angles[3]:
                    angle_result = validate_angles(
                        fixed_chains, moving_chains, c4_axes[1:5], [angle1, angle2, angle3, angle4]
                    )
                    if angle_result:
                        valid_combinations.append(angle_result)
    if valid_combinations:
        valid_combinations.sort(key=lambda x: x["rmsd"])
        all_chains = valid_combinations[0]["chains"]
        fixed_chains = [
            all_chains[6], all_chains[17], all_chains[12], all_chains[11]
        ]
        moving_chains = new_chains[20: 24]
        final_valid_angles = []
        for step in range(num_steps):
            angle_result = check_rotated_angle_2(
                fixed_chains, moving_chains, c4_axes[5], step * angle_step
            )
            if angle_result:
                final_valid_angles.append(angle_result)
        if final_valid_angles:
            final_valid_angles.sort(key=lambda x: x["rmsd"])
            all_chains.extend(final_valid_angles[0]["chains"])
            initial_result = {
                "c4_idx": c4_idx, "chains": all_chains, "c4_axes": c4_axes,
                "sliding_range": sliding_range, "clash_count": clash_count
            }
            return initial_result
        else:
            return None
    else:
        return None

# traverse all C4 symmetric modes and retain those can form O-symmetry
def generate_initial_results_2(c4_modes_all, num_cpus, angle_step):
    initial_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(copy_c4_chains, i+1, c4_mode, angle_step)
            for i, c4_mode in enumerate(c4_modes_all)
        ]
        for future in as_completed(futures):
            result = future.result()
            if result:
                initial_results.append(result)
    gc.collect()

    if initial_results:
        change_chain_ids(initial_results)
    return initial_results

# perform sliding optimization on the C4-formed O-symmetric systems
def sliding_optimization_2(idx, result, md_steps):
    c4_idx = result["c4_idx"]
    c4_axes = result["c4_axes"]
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
        for i, c4_axis in enumerate(c4_axes):
            for chain in start_chains[4*i: 4*(i+1)]:
                new_chain = chain.copy()
                for atom in new_chain.get_atoms():
                    atom.coord += sliding_step * c4_axis
                new_chains.append(new_chain)
        total_clash_count = count_clashes(new_chains)
        sliding_results.append(
            {"c2_idx": 0, "c3_idx": 0, "c4_idx": c4_idx,
             "step": step, "chains": new_chains, "clash": total_clash_count}
        )
    contact_indices = []
    for i, sliding_result in enumerate(sliding_results):
        if sliding_result["clash"] > 6 * clash_count:
            contact_indices.append(i)
    sliding_results = sliding_results[contact_indices[0]-1: contact_indices[-1]+2]
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

# optimize all C4-formed O-symmetric systems
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


############    MAIN    ############
def main():
    parser = ArgumentParser()
    parser.add_argument("-n2", "--n2", type=int, default=10)
    parser.add_argument("-n3", "--n3", type=int, default=10)
    parser.add_argument("-n4", "--n4", type=int, default=10)
    parser.add_argument("-w", "--workers", type=int, default=20)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("-m", "--md", type=int, default=500)
    parser.add_argument("-n", "--nmax", type=int, default=10)

    args = parser.parse_args()
    num_c2 = args.n2
    num_c3 = args.n3
    num_c4 = args.n4
    max_workers = args.workers
    num_cpus = get_num_cpus(max_workers)
    output_file = args.output
    o_stem, o_suffix = output_file.rsplit(".", 1)
    md_steps = args.md
    num_models = args.nmax

    start_time = time.time()
    # print(f"[INFO] Start O symmetric docking", flush=True)

    # Local symmetric modes sampling
    parent_dir = os.path.dirname(current_dir)
    mono_file = "A.pdb"
    cn_out_files = ["A_c2.out", "A_c4.out"]
    cn_pdb_files = ["A2.pdb", "A4.pdb"]
    cn_log_files = ["A2.log", "A4.log"]
    cn_types = [2, 4]
    cn_nums = [num_c2, num_c4]
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

    stem = "A2"
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(run_csym_dock, stem, i, 3, num_c3) for i in range(1, 11)
        ]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"[ERROR] Failed in local symmetry sampling: {e}")
                sys.exit(1)
    gc.collect()

    # Octahedral symmetry construction
    c23_modes_all = []
    c4_modes_all = []
    for i in range(1, num_c2+1):
        for j in range(1, num_c3+1):
            pdb_file = f"A23_{i}_{j}.pdb"
            c23_mode = load_c23_symmetry_mode(pdb_file)
            c23_modes_all.append(
                {"c2_idx": i, "c3_idx": j, "chains": c23_mode}
            )
    for i in range(1, num_c4+1):
        pdb_file = f"A4_{i}.pdb"
        c4_mode = load_c4_symmetry_mode(pdb_file)
        c4_modes_all.append(c4_mode)

    initial_results = generate_initial_results_1(c23_modes_all, c4_modes_all, num_cpus)
    if initial_results:
        selected_results = []
        c23_results = {}
        for result in initial_results:
            c2_idx, c3_idx = result["c2_idx"], result["c3_idx"]
            c23_indices = (c2_idx, c3_idx)
            if c23_indices not in c23_results:
                c23_results[c23_indices] = []
            c23_results[c23_indices].append(result)
        for c23_indices, results in c23_results.items():
            if not results:
                continue
            results.sort(key=lambda x: x["c4_idx"])
            selected_results.append(results[0])
        refined_results = refine_results_1(selected_results, num_cpus, md_steps, num_models)
        end_time = time.time()
        with open(f"{current_dir}/aa2.log", 'w+') as fout:
            for i, result in enumerate(refined_results):
                j, k, m = result["c2_idx"], result["c3_idx"], result["c4_idx"]
                it_score = result["it_score"]
                print(f"Model {i + 1}, ({j}, {k})-{m}, IT-score: {it_score}", file=fout)
                model_file = f"{current_dir}/{o_stem}_{i + 1}.{o_suffix}"
                shutil.copy(result["file"], model_file)
            print(f"\nOctahedral symmetric docking cost {end_time - start_time}",
                  file=fout, flush=True)

    else:
        print("No valid C2-C3-C4 combination, try to construct Osymmetry manually.")
        angle_step = np.pi / 12
        new_initial_results = generate_initial_results_2(c4_modes_all, num_cpus, angle_step)
        if new_initial_results:
            refined_results = refine_results_2(new_initial_results, num_cpus, md_steps, num_models)
            end_time = time.time()
            with open(f"{current_dir}/aa2.log", 'w+') as fout:
                for i, result in enumerate(refined_results):
                    j, k, m = result["c2_idx"], result["c3_idx"], result["c4_idx"]
                    it_score = result["it_score"]
                    print(f"Model {i + 1}, ({j}, {k})-{m}, IT-score: {it_score}", file=fout)
                    model_file = f"{current_dir}/{o_stem}_{i + 1}.{o_suffix}"
                    shutil.copy(result["file"], model_file)
                print(f"\nOctahedral symmetric docking cost {end_time - start_time}",
                      file=fout, flush=True)
            print_results(refined_results, output_file, num_cpus, md_steps)
        else:
            print(f"[WARN] No valid Osymmetry configurations found.")
            print(f"[WARN] Try more C4 tetramers or smaller angle steps.")
            end_time = time.time()

    shutil.rmtree(temp_dir)
    print(f"\nOctahedral symmetric docking cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    main()
