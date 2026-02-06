import os
import gc
import json
import sys
from pathlib import Path
import subprocess
from parse_pdb import PDBParser, load_format_line, write_pdb_file
from scipy.spatial.transform import Rotation
import copy
import numpy as np
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil
import time

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
pdb_parser = PDBParser()
format_lines = load_format_line(os.path.join(BASE_PATH, "format.pdb"))
hdock = os.path.join(TOOL_PATH, "hdock")
chdock = os.path.join(TOOL_PATH, "chdock")
compcn = os.path.join(TOOL_PATH, "compcn")
splitmodels_script = os.path.join(TOOL_PATH, "splitmodels.py")

steric_clash_dist = 2


# select appropriate number of parallel CPU cores
def get_num_cpus(max_workers: int):
    cpu_cap = int(os.getenv("NUM_CPUS", max_workers))
    cpu_count = os.cpu_count()
    return max(1, min(cpu_count, cpu_cap))

# perform factorization on the integer n
def factorize(n):
    factors = []
    for i in range(2, n + 1):
        if n % i == 0:
            factors.append(i)
    return factors

# load hdock output file
def load_out_information(out_file, num_lines: int = 5):
    with open(out_file, 'r') as file:
        lines = file.readlines()
    if len(lines) < 5 + num_lines:
        raise ValueError(f"[ERROR] Check input hdock out file: {out_file}")
    line = lines[2]
    parts = line.strip().split()
    rands = [float(parts[i]) for i in range(2, 5)]
    line = lines[3]
    parts = line.strip().split()
    rec_file = parts[0]
    line = lines[4]
    parts = line.strip().split()
    lig_file = parts[0]
    lig_centroid = np.array([float(parts[i]) for i in range(1, 4)])
    solutions = []
    for i in range(5, 5 + num_lines):
        line = lines[i]
        parts = line.strip().split()
        solution = [float(parts[j]) for j in range(6)]
        solutions.append(solution)
    return rands, rec_file, lig_file, lig_centroid, solutions

# calculate rotation matrix from Euler angles
def calculate_rot(angles):
    phi, theta, psi = angles[:3]
    cpsi = math.cos(psi)
    spsi = math.sin(psi)
    cthe = math.cos(theta)
    sthe = math.sin(theta)
    cphi = math.cos(phi)
    sphi = math.sin(phi)

    rot = [[0.0] * 3 for _ in range(3)]
    rot[0][0] = cpsi * cphi - cthe * sphi * spsi
    rot[0][1] = cpsi * sphi + cthe * cphi * spsi
    rot[0][2] = spsi * sthe
    rot[1][0] = -spsi * cphi - cthe * sphi * cpsi
    rot[1][1] = -spsi * sphi + cthe * cphi * cpsi
    rot[1][2] = cpsi * sthe
    rot[2][0] = sthe * sphi
    rot[2][1] = -sthe * cphi
    rot[2][2] = cthe
    return np.array(rot)

# apply transformation to the chain based on the given angles and translations
def transform_chain(chain, centroid, rands, solution):
    new_chain = chain.copy()
    rot = calculate_rot(rands)
    for atom in new_chain.get_atoms():
        atom.coord = rot @ (atom.coord - centroid)
    angles = solution[:3]
    rot = calculate_rot(angles)
    trans = np.array(solution[3:])
    for atom in new_chain.get_atoms():
        atom.coord = rot @ atom.coord + trans + centroid
    return new_chain

# calculate rotation matrix for aligning two sets of coordinates using Kabsch's algorithm
def kabsch_rot(coords1, coords2):
    H = np.dot(coords1.T, coords2)
    U, S, VT = np.linalg.svd(H)
    rot = np.dot(VT.T, U.T)
    if np.linalg.det(rot) < 0.0:
        VT[2, :] *= -1
        rot = np.dot(VT.T, U.T)
    return rot

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
        raise ValueError("[ERROR] The sizes of input chains are different")

    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    centered_coords1 = coords1 - centroid1
    centered_coords2 = coords2 - centroid2
    rot = kabsch_rot(centered_coords2, centered_coords1)
    trans = centroid2 - np.dot(centroid1, rot)
    return rot, trans

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

# calculate CA rmsd between two chains
def calculate_ca_rmsd(chain1, chain2):
    coords1 = np.array([atom.coord for atom in chain1.get_atoms() if atom.name == "CA"])
    coords2 = np.array([atom.coord for atom in chain2.get_atoms() if atom.name == "CA"])
    diff = coords1 - coords2
    ca_rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return ca_rmsd

# calculate number of spatial conflicts in the chain group
def count_clashes(chains):
    chains_coords = [np.array([atom.coord for atom in chain.get_atoms()]) for chain in chains]
    total_clash_count = 0
    for i in range(len(chains_coords)):
        for j in range(i + 1, len(chains_coords)):
            coords1, coords2 = chains_coords[i], chains_coords[j]
            distances = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=2)
            clash_count = np.sum(distances < steric_clash_dist)
            total_clash_count += clash_count
    return total_clash_count

# check whether the two given chains satisfy C2 symmetry
def check_c2_symmetry_rmsd(chain1, chain2, rmsd_threshold):
    coords1 = [atom.coord for atom in chain1.get_atoms()]
    coords2 = [atom.coord for atom in chain2.get_atoms()]
    if len(coords1) != len(coords2):
        raise ValueError("[ERROR] Length of input chains is different")

    coords = coords1 + coords2
    centroid = np.mean(coords, axis=0)
    atm_idx = len(coords1) // 2
    c2_axis = np.cross(coords1[0] - coords2[0], coords1[atm_idx] - coords2[atm_idx])
    c2_axis /= np.linalg.norm(c2_axis)
    angle = np.pi
    rotation_matrix = Rotation.from_rotvec(angle * c2_axis).as_matrix()
    rotated_chain = chain2.copy()
    for atom in rotated_chain.get_atoms():
        atom.coord = rotation_matrix @ (atom.coord - centroid) + centroid
    ca_rmsd = calculate_ca_rmsd(chain1, rotated_chain)
    return ca_rmsd < rmsd_threshold

# check whether the transformation extracted from the binary chain group
# can form cyclic symmetry
def search_csym_for_chains(idx, start_chains, chain_num, tolerance, rmsd_threshold):
    if len(start_chains) != 2:
        raise ValueError("[ERROR] Check input chains num")
    c2_symmetry = check_c2_symmetry_rmsd(start_chains[0], start_chains[1], rmsd_threshold)
    if c2_symmetry:
        return {"idx": idx, "cn_symmetry": 2}

    rot, trans = calculate_align_transformation(start_chains[0], start_chains[1])
    clash_count = count_clashes(start_chains)
    cn_symmetry = False
    for n in range(3, chain_num + 1):
        chains = copy.deepcopy(start_chains)
        for i in range(n - 2):
            new_chain = chains[-1].copy()
            for atom in new_chain.get_atoms():
                atom.coord = np.dot(atom.coord, rot) + trans
            chains.append(new_chain)
        total_clash_count = count_clashes(chains)
        if total_clash_count > (n - 1) * clash_count + tolerance:
            break
        chains = center_chains(chains)
        centroids = [calculate_centroid(chain) for chain in chains]
        cn_axis = np.cross(centroids[-1] - centroids[0], centroids[1] - centroids[0])
        cn_axis = cn_axis / np.linalg.norm(cn_axis)
        angle = 2 * np.pi / n
        rotation_matrix = Rotation.from_rotvec(angle * cn_axis).as_matrix()
        all_chains_valid = True
        for i, chain in enumerate(chains):
            chain1 = chains[(i - 1) % n]
            chain2 = chains[(i + 1) % n]
            rotated_chain = chain.copy()
            for atom in rotated_chain.get_atoms():
                atom.coord = rotation_matrix @ atom.coord
            ca_rmsd = min(
                calculate_ca_rmsd(chain1, rotated_chain), calculate_ca_rmsd(chain2, rotated_chain)
            )
            if ca_rmsd >= rmsd_threshold:
                all_chains_valid = False
                break
        if all_chains_valid:
            cn_symmetry = True
            break
    return {"idx": idx, "cn_symmetry": n} if cn_symmetry else None

# compare the CHDOCK sampled Cn symmetric structure with
# the Cn symmetric structure derived from the pairwise docking transformation
def match_cn_symmetry(csym_chains, cn_chains, rmsd_threshold):
    cn_idx = cn_chains["cn_idx"]
    temp_chains = cn_chains["chains"]
    temp_file = cn_chains["file"]
    rot, trans = calculate_align_transformation(csym_chains[0], temp_chains[0])
    new_chain = csym_chains[1].copy()
    for atom in new_chain.get_atoms():
        atom.transform(rot, trans)
    if len(temp_chains) == 2:
        ca_rmsd = calculate_ca_rmsd(new_chain, temp_chains[1])
    else:
        ca_rmsd = min(
            calculate_ca_rmsd(new_chain, temp_chains[-1]),
            calculate_ca_rmsd(new_chain, temp_chains[1])
        )
    if ca_rmsd < rmsd_threshold:
        match_result = {"cn_idx": cn_idx, "file": temp_file}
        return match_result
    else:
        return None

# perform pairwise docking on monomer using HDOCK;
# and check whether the top 10 transformations from pairwise docking can form cyclic symmetry;
# use CHDOCK to further verify the cyclic-like transformation from monomer pairwise docking
def check_csym_from_monodock(chain_num, current_path, num_cpus, rmsd_threshold):
    # pairwise docking
    os.chdir(current_path)
    mono_file = "A.pdb"
    mono_out_file = "AA.out"
    subprocess.run(f"{hdock} {mono_file} {mono_file} -out {mono_out_file} > /dev/null",
                   shell=True, check=True)
    rands, rec_file, lig_file, lig_centroid, solutions = load_out_information(mono_out_file, 10)
    if not (Path(rec_file).exists() and Path(lig_file).exists()):
        raise FileNotFoundError(f"[ERROR] Please check path for rec: {rec_file} or lig: {lig_file}")

    chain1 = pdb_parser.get_structure("rec", rec_file)[0]["A"]
    chain2 = pdb_parser.get_structure("lig", lig_file)[0]["A"]
    all_start_chains = []
    for solution in solutions:
        start_chains = [chain1]
        new_chain2 = transform_chain(chain2, lig_centroid, rands, solution)
        new_chain2.id = "B"
        start_chains.append(new_chain2)
        all_start_chains.append(start_chains)

    csym_results = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(
                search_csym_for_chains, i, start_chains, chain_num, 50, rmsd_threshold
            ) for i, start_chains in enumerate(all_start_chains)
        ]
        for future in as_completed(futures):
            result = future.result()
            if result:
                csym_results.append(result)

    # further validation by Cn construction using CHDOCK
    if csym_results:
        csym_results.sort(key=lambda x: x["cn_symmetry"], reverse=True)
        out_idx = csym_results[0]["idx"]
        cn_symmetry = csym_results[0]["cn_symmetry"]
        if chain_num % cn_symmetry == 0:
            mono_csym_chains = all_start_chains[out_idx]
            cn_out_file = f"A_c{cn_symmetry}.out"
            cn_pdb_file = f"A{cn_symmetry}.pdb"
            cn_log_file = f"A{cn_symmetry}.log"
            subprocess.run(
                f"{chdock} {mono_file} {mono_file} -symmetry {cn_symmetry} -out {cn_out_file} > {cn_log_file} &&"
                f"{compcn} {cn_out_file} {cn_pdb_file} -nmax 10 -complex > /dev/null &&"
                f"python {splitmodels_script} {cn_pdb_file}", shell=True, check=True
            )
            cn_chains_all = []
            for i in range(1, 11):
                pdb_file = f"A{cn_symmetry}_{i}.pdb"
                model = pdb_parser.get_structure("", pdb_file)[0]
                chains = list(model.get_chains())
                cn_chains_all.append(
                    {"cn_idx": i, "chains": chains, "file": pdb_file}
                )
            match_results_all = []
            with ProcessPoolExecutor(max_workers=num_cpus) as executor:
                futures = [
                    executor.submit(
                        match_cn_symmetry, mono_csym_chains, cn_chains, rmsd_threshold
                    ) for cn_chains in cn_chains_all
                ]
                for future in as_completed(futures):
                    result = future.result()
                    if result:
                        match_results_all.append(result)
            gc.collect()

            if match_results_all:
                match_results_all.sort(key=lambda x: x["cn_idx"])
                selected_cn_result = match_results_all[0]
                selected_cn_pdb_file = selected_cn_result["file"]
                return cn_symmetry, selected_cn_pdb_file, cn_out_file
            else:
                return 1, mono_file, mono_out_file
        else:
            return 1, mono_file, mono_out_file
    else:
        return 1, mono_file, mono_out_file


####################      MAIN      ####################
def main():
    file_stoi = sys.argv[1]

    current_dir = os.path.dirname(os.path.abspath(file_stoi))
    # num_cpus = int(os.getenv("NUM_CPUS", 20))
    num_cpus = get_num_cpus(20)

    start_time = time.time()

    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    build_types = ["asym", "assemble"]
    if len(stoi_data) > 1:
        if len(stoi_data) >= 5:
            build_types = ["assemble"]
        print("\nHeteromeric complex")
        with open(f"{current_dir}/build_types.txt", 'w+') as fout:
            for i in build_types:
                print(i, file=fout)
                model_path = f"{current_dir}/{i}/"
                Path(model_path).mkdir(exist_ok=True, parents=True)
                if i == "asym":
                    unique_chain_ids = [j["asym_ids"][0] for j in stoi_data]
                    for chain_id in unique_chain_ids:
                        shutil.copy(f"{current_dir}/{chain_id}.pdb", model_path)
    else:
        print("\nHomomeric complex")
        rmsd_threshold = 2.0
        chain_num = len(stoi_data[0]["asym_ids"])
        if chain_num == 12:
            build_types.append("cubicT")
        elif chain_num == 24:
            build_types.append("cubicO")
        cn, cn_pdb_file, cn_out_file = check_csym_from_monodock(
            chain_num, current_dir, num_cpus, rmsd_threshold
        )
        if cn == 1:
            factors = factorize(chain_num)
            for i in factors:
                build_types.append(f"c{i}")
                if i == 2:
                    j = chain_num // 2
                    build_types.append(f"d{j}")

            with open(f"{current_dir}/build_types.txt", 'w+') as fout:
                for i in build_types:
                    print(i, file=fout)
                    model_path = f"{current_dir}/{i}/"
                    Path(model_path).mkdir(exist_ok=True, parents=True)
                    if i == "asym":
                        shutil.copy(cn_pdb_file, model_path)
                        output_path = f"{model_path}/pairwise_out/"
                        Path(output_path).mkdir(exist_ok=True, parents=True)
                        shutil.copy(cn_out_file, output_path)
                    elif i == "assemble":
                        continue
                    else:
                        shutil.copy(cn_pdb_file, model_path)

        else:
            if cn == chain_num:
                build_types.append(f"c{cn}")
            else:
                complementary_cn = chain_num // cn
                if cn == 2 or complementary_cn == 2:
                    build_types.extend(
                        [f"c{cn}", f"c{complementary_cn}", f"d{max(cn, complementary_cn)}"]
                    )
                else:
                    build_types.extend([f"c{cn}", f"c{complementary_cn}"])

            with open(f"{current_dir}/build_types.txt", 'w+') as fout:
                for i in build_types:
                    print(i, file=fout)
                    model_path = f"{current_dir}/{i}/"
                    Path(model_path).mkdir(exist_ok=True, parents=True)
                    if i == "asym":
                        shutil.copy(cn_pdb_file, model_path)
                    elif i == "assemble":
                        continue
                    else:
                        shutil.copy(f"{current_dir}/A.pdb", model_path)

    end_time = time.time()
    print(f"Determining modeling strategies cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    main()
