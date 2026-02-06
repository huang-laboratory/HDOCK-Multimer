import sys
import json
import os
from pathlib import Path
from parse_pdb import PDBParser, Chain
import numpy as np
import shutil
import time


protein_letters = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}
pdb_parser = PDBParser()
all_chain_ids = []
homo_chains_list = []


############    Utility Functions    ############
# collect models from multi-body docking (asymmetric docking)
def collect_multi_dock_models(current_dir):
    report_log = f"{current_dir}/asym/multi.out"
    if Path(report_log).is_file():
        with open(report_log, 'r') as file:
            lines = file.readlines()
            parts = lines[-1].strip().split()
            it_score = float(parts[-1])
        model_file = f"{current_dir}/asym/models/multi.pdb"
        multi_dock_result = {
            "it_score": it_score, "file": model_file, "type": "multi_dock"
        }
        return multi_dock_result
    else:
        raise FileNotFoundError(f"[ERROR] Please check {current_dir}/asym/")

# collect models from assembly
def collect_assemble_models(current_dir):
    report_log = f"{current_dir}/assemble/aa2.log"
    assemble_results = []
    assemblied_all = False
    if Path(report_log).is_file():
        with open(report_log, 'r') as file:
            for line in file:
                if line.startswith("Iteration"):
                    assemblied_all = True
                    parts = line.strip().split()
                    it_score = float(parts[-1])
                    iteration, num = parts[1], parts[3]
                    model_file = f"{current_dir}/assemble/assemble{iteration}_{num}.pdb"
                    if it_score != 0.0:
                        assemble_results.append(
                            {"it_score": it_score, "file": model_file, "flag": "assemble"}
                        )
                    else:
                        print(f"[WARN] Please check {model_file}", flush=True)
    else:
        raise FileNotFoundError(f"[ERROR] Please check {current_dir}/assemble/")

    if assemblied_all:
        return assemble_results, assemblied_all
    else:
        report_log_docking = f"{current_dir}/assemble/dock/multi.out"
        if Path(report_log_docking).is_file():
            with open(report_log_docking, 'r') as file:
                lines = file.readlines()
                parts = lines[-1].strip().split()
                it_score = float(parts[-1])
                model_file = f"{current_dir}/assemble/dock/models/multi.pdb"
                assemble_results.append(
                    {"it_score": it_score, "file": model_file, "flag": "assemble_dock"}
                )
            return assemble_results, assemblied_all
        else:
            raise FileNotFoundError(f"[ERROR] Please check {current_dir}/assemble/dock/")

# collect models from T-symmetric / O-symmetric docking
def collect_cubic_models(current_dir, build_type):
    cubic_type = build_type[-1]
    cubic_results = []
    report_log = f"{current_dir}/{build_type}/aa2.log"
    if Path(report_log).is_file():
        with open(report_log, 'r') as file:
            for line in file:
                if line.startswith("Model"):
                    parts = line.strip().split()
                    num, it_score = parts[1][:-1], float(parts[-1])
                    model_file = f"{current_dir}/{build_type}/{cubic_type}sym_{num}.pdb"
                    if it_score != 0.0:
                        cubic_results.append(
                            {"it_score": it_score, "file": model_file, "flag": build_type}
                        )
                    else:
                        print(f"[WARN] Please check {model_file}", flush=True)
        return cubic_results
    else:
        raise FileNotFoundError(f"[ERROR] Please check {current_dir}/{build_type}/")

# collect models from cyclic symmetric docking
def collect_cn_models(current_dir, build_type):
    report_log = f"{current_dir}/{build_type}/aa2.log"
    cn_results = []
    if Path(report_log).is_file():
        with open(report_log, 'r') as file:
            for line in file:
                if line.startswith("Model"):
                    parts = line.strip().split()
                    num, it_score = parts[1][:-1], float(parts[-1])
                    model_file = f"{current_dir}/{build_type}/cn_{num}.pdb"
                    cn_results.append(
                        {"it_score": it_score, "file": model_file, "flag": build_type}
                    )
        return cn_results
    else:
        raise FileNotFoundError(f"[ERROR] Please check {current_dir}/{build_type}/")

# collect models from dihedral symmetric docking
def collect_dn_models(current_dir, build_type):
    report_log = f"{current_dir}/{build_type}/aa2.log"
    dn_results = []
    if Path(report_log).is_file():
        with open(report_log, 'r') as file:
            for line in file:
                if line.startswith("Model"):
                    parts = line.strip().split()
                    num, it_score = parts[1][:-1], float(parts[-1])
                    model_file = f"{current_dir}/{build_type}/dn_{num}.pdb"
                    dn_results.append(
                        {"it_score": it_score, "file": model_file, "flag": build_type}
                    )
        return dn_results
    else:
        raise FileNotFoundError(f"[ERROR] Please check {current_dir}/{build_type}/")

# determine the range of residues for calculating RMSD
def determine_res_indices(current_dir, stoi_data, assemble_flag, dock_flag):
    global all_chain_ids, homo_chains_list
    res_indices_mono = {}
    res_indices_assemble = {}
    final_res_indices = {}
    unique_chain_ids = {}
    all_chain_ids = []
    homo_chains_list = []
    for i in stoi_data:
        asym_ids, sequence = i["asym_ids"], i["sequence"]
        if len(asym_ids) > 1:
            homo_chains_list.append(asym_ids)
        all_chain_ids.extend(asym_ids)
        unique_chain_ids[asym_ids[0]] = asym_ids
        for chain_id in asym_ids:
            res_indices_assemble[chain_id] = (1, len(sequence))
    if not assemble_flag:
        file_res_indices1 = f"{current_dir}/assemble/res_indices_A.txt"
        if Path(file_res_indices1).is_file():
            with open(file_res_indices1, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    res_indices_assemble[parts[0]] = (int(parts[1]), int(parts[2]))
        else:
            raise FileNotFoundError(f"[ERROR] Missing file: {file_res_indices1}")

    if dock_flag:
        file_res_indices2 = f"{current_dir}/res_indices_D.txt"
        unique_res_indices = {}
        if Path(file_res_indices2).exists():
            with open(file_res_indices2, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    unique_res_indices[parts[0]] = (int(parts[1]), int(parts[2]))
        else:
            raise FileNotFoundError(f"[ERROR] Missing file: {file_res_indices2}")
        for chain_id, asym_ids in unique_chain_ids.items():
            for asym_id in asym_ids:
                res_indices_mono[asym_id] = unique_res_indices[chain_id]

        all_chain_ids.sort()
        for chain_id in all_chain_ids:
            start_idx1, end_idx1 = res_indices_assemble[chain_id]
            start_idx2, end_idx2 = res_indices_mono[chain_id]
            start_idx = max(start_idx1, start_idx2) + 1
            end_idx = min(end_idx1, end_idx2) - 1
            if start_idx > end_idx:
                raise ValueError(f"[ERROR] Wrong residue indices for chain {current_dir}/{chain_id}.pdb")
            final_res_indices[chain_id] = (start_idx, end_idx)
    else:
        final_res_indices = res_indices_assemble
    return final_res_indices

# filter the chain based on the given residue index range
def filter_chain(chain, start_idx, end_idx):
    filtered_chain = Chain(chain.id)
    for residue in chain.get_residues():
        if residue.resseq in range(start_idx, end_idx + 1):
            filtered_chain.residues.append(residue)
    return filtered_chain

# calculate centroid for given coordinates
def calculate_centroid(coords):
    return np.mean(coords, axis=0)

# calculate centroid for given chains
def calculate_centroid_chains(chains):
    if isinstance(chains, list):
        coords = [atom.coord for chain in chains for atom in chain.get_atoms()]
    else:
        coords = [atom.coord for atom in chains.get_atoms()]
    centroid = calculate_centroid(coords)
    return centroid

# calculate rotation matrix for aligning two sets of coordinates using Kabsch's algorithm
def kabsch_rot(coords1, coords2):
    H = np.dot(coords1.T, coords2)
    U, S, VT = np.linalg.svd(H)
    rot = np.dot(VT.T, U.T)
    if np.linalg.det(rot) < 0.0:
        VT[2, :] *= -1
        rot = np.dot(VT.T, U.T)
    return rot

# calculate CA rmsd for two sets of coordinates
def calculate_rmsd_coords(coords1, coords2):
    centroid1 = calculate_centroid(coords1)
    centroid2 = calculate_centroid(coords2)
    coords1 -= centroid1
    coords2 -= centroid2
    rot = kabsch_rot(coords1, coords2)
    coords2 = np.dot(coords2, rot)
    diff = coords1 - coords2
    ca_rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return ca_rmsd

# calculate CA rmsd for two complexes (optimal subunits correspondence)
def calculate_rmsd_orientation(orientation1, orientation2, file1, file2):
    coords1 = np.array([
        residue["CA"].get_coord() for chain in orientation1.values()
        for residue in chain.get_residues() if "CA" in residue.atoms
    ])
    coords2 = np.array([
        residue["CA"].get_coord() for chain in orientation2.values()
        for residue in chain.get_residues() if "CA" in residue.atoms
    ])
    if coords1.shape != coords2.shape:
        raise ValueError(f"[ERROR] Not equal CA atoms number in {file1}, {file2}")
    ca_rmsd = calculate_rmsd_coords(coords1, coords2)
    return ca_rmsd

# search for the optimal correspondence between subunits from two complexes
def search_subunits_mapping(orientation1, orientation2):
    centroids1 = {chain_id: calculate_centroid_chains(chain) for chain_id, chain in orientation1.items()}
    centroids2 = {chain_id: calculate_centroid_chains(chain) for chain_id, chain in orientation2.items()}

    best_mapping = {chain_id: chain_id for chain_id in all_chain_ids}
    max_iter = 10
    for _ in range(max_iter):
        improved = False
        for homo_chains in homo_chains_list:
            homo_best_mapping = {chain_id: best_mapping[chain_id] for chain_id in homo_chains}
            homo_centroids1 = np.array([centroids1[chain_id] for chain_id in homo_chains])
            homo_centroids2 = np.array([
                centroids2[homo_best_mapping[chain_id]] for chain_id in homo_chains
            ])
            best_rmsd = calculate_rmsd_coords(homo_centroids1, homo_centroids2)
            for i, chain1 in enumerate(homo_chains):
                for j, chain2 in enumerate(homo_chains):
                    if i >= j:
                        continue
                    swapped_mapping = homo_best_mapping.copy()
                    swapped_mapping[chain1], swapped_mapping[chain2] = swapped_mapping[chain2], swapped_mapping[chain1]
                    homo_centroids2 = np.array([
                        centroids2[swapped_mapping[chain_id]] for chain_id in homo_chains
                    ])
                    current_rmsd = calculate_rmsd_coords(homo_centroids1, homo_centroids2)
                    if current_rmsd < best_rmsd:
                        best_rmsd = current_rmsd
                        homo_best_mapping = swapped_mapping
                        improved = True
            if improved:
                for chain_id in homo_chains:
                    best_mapping[chain_id] = homo_best_mapping[chain_id]
        if not improved:
            break

    updated_orientation2 = {
        chain_id: orientation2[best_mapping[chain_id]].copy() for chain_id in all_chain_ids
    }
    for chain_id in updated_orientation2:
        updated_orientation2[chain_id].id = chain_id
    return updated_orientation2

# calculate CA rmsd for two complexes
def calculate_rmsd(orientation1, orientation2, file1, file2):
    if homo_chains_list:
        orientation2 = search_subunits_mapping(orientation1, orientation2)
        ca_rmsd = calculate_rmsd_orientation(orientation1, orientation2, file1, file2)
    else:
        ca_rmsd = calculate_rmsd_orientation(orientation1, orientation2, file1, file2)
    return ca_rmsd


############    MAIN    ############
def main():
    file_stoi = sys.argv[1]
    nmax = int(sys.argv[2])

    start_time = time.time()
    current_dir = os.path.dirname(os.path.abspath(file_stoi))
    print(current_dir)

    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    # read model files of different modeling strategies
    build_types_file = f"{current_dir}/build_types.txt"
    build_types_list = []
    if Path(build_types_file).is_file():
        with open(build_types_file, 'r') as file:
            for line in file:
                build_type = line.strip()
                build_types_list.append(build_type)
    else:
        raise FileNotFoundError(f"[ERROR] Missing build_type.txt, please check {current_dir}/")

    all_results = []
    assemblied_all = True
    dock_performed = False
    for build_type in build_types_list:
        if build_type == "assemble":
            assemble_results, assemblied_all = collect_assemble_models(current_dir)
            all_results.extend(assemble_results)
        elif build_type == "asym":
            multi_dock_result = collect_multi_dock_models(current_dir)
            all_results.append(multi_dock_result)
            dock_performed = True
        elif build_type.startswith("cubic"):
            cubic_results = collect_cubic_models(current_dir, build_type)
            all_results.extend(cubic_results)
        elif build_type.startswith("c"):
            cn_results = collect_cn_models(current_dir, build_type)
            all_results.extend(cn_results)
        elif build_type.startswith("d"):
            dn_results = collect_dn_models(current_dir, build_type)
            all_results.extend(dn_results)
        else:
            raise ValueError(f"[ERROR] Unsupported building strategy: {build_type}")

    # models clustering
    res_indices = determine_res_indices(current_dir, stoi_data, assemblied_all, dock_performed)
    for result in all_results:
        model_file = result["file"]
        model = pdb_parser.get_structure("", model_file)[0]
        orientation = {}
        for chain_id in all_chain_ids:
            start_idx, end_idx = res_indices[chain_id]
            chain = model[chain_id]
            filtered_chain = filter_chain(chain, start_idx, end_idx)
            orientation[chain_id] = filtered_chain
        result["orientation"] = orientation

    all_results.sort(key=lambda x: x["it_score"])
    clstr_results = []
    rmsd_threshold = 5.0
    for result1 in all_results:
        added_to_clstr = False
        for result2 in clstr_results:
            ca_rmsd = calculate_rmsd(
                result1["orientation"], result2["orientation"], result1["file"], result2["file"]
            )
            if ca_rmsd < rmsd_threshold:
                added_to_clstr = True
                break
        if not added_to_clstr:
            clstr_results.append(result1)
        if len(clstr_results) >= nmax:
            break

    end_time = time.time()
    result_dir = f"{current_dir}/results/"
    Path(result_dir).mkdir(exist_ok=True, parents=True)
    with open(f"{result_dir}/results.log", 'w+') as fout:
        for i, result in enumerate(clstr_results):
            file1 = result["file"]
            file2 = f"{result_dir}/model_{i + 1}.pdb"
            it_score = result["it_score"]
            shutil.copy(file1, file2)
            print(f"Model {i + 1} IT-score: {it_score}", file=fout, flush=True)
        print(f"\nModels clustering cost {end_time - start_time}", file=fout, flush=True)

    print(f"\nModels clustering cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    main()
