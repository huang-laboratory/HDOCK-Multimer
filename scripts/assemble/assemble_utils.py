# define tool functions used in the assembly strategy

import numpy as np
import pandas as pd
from pathlib import Path
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc
import time
import config
from config import managed_tempfile, current_dir
import sys
import shutil

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
openmm_script = os.path.join(TOOL_PATH, "refine_model_openmm.py")
scorecom = os.path.join(TOOL_PATH, "scorecom.sh")
sys.path.append(BASE_PATH)
from parse_pdb import Chain, write_pdb_file


############    Utility Functions for Assembly    ############
# generate inverse for a given transformation (rotation + translation)
def reverse_transform(rot, trans):
    inverse_rot = rot.T
    inverse_trans = -np.dot(trans, inverse_rot)
    return inverse_rot, inverse_trans

# generate inverse for a list of transformations
def reverse_transform_list(old_transforms: list):
    inverse_transforms = []
    for rot, trans in reversed(old_transforms):
        inverse_rot, inverse_trans = reverse_transform(rot, trans)
        inverse_transforms.append((inverse_rot, inverse_trans))
    return inverse_transforms

# compress a list of transformations into one
def combine_transform_list(transforms: list):
    total_rot = np.eye(3)
    total_trans = np.zeros(3)
    for rot, trans in transforms:
        total_rot = np.dot(total_rot, rot)
        total_trans = np.dot(total_trans, rot) + trans
    return total_rot, total_trans

# apply a list of transformation to the chain
def apply_transform_list(chain: Chain, transforms: list):
    total_rot, total_trans = combine_transform_list(transforms)
    for atom in chain.get_atoms():
        atom.transform(total_rot, total_trans)

# calculate IT-score for the complex directly
def calculate_itscore(orientation, i):
    with managed_tempfile(".pdb", delete_ok=False) as model_file:
        write_pdb_file(orientation, model_file, config.format_lines)
        result = subprocess.run(f"bash {scorecom} {model_file}", shell=True, check=True,
                                stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        it_score = float(result.stdout.strip())
        return i, it_score, model_file

# refine structure using OpenMM energy minimization, then calculate IT-score
def calculate_itscore_md(orientation, i, md_steps):
    with managed_tempfile(".pdb") as model_file, \
            managed_tempfile(".pdb", delete_ok=False) as md_file:
        write_pdb_file(orientation, model_file, config.format_lines)
        subprocess.run(f"python {openmm_script} {model_file} {md_steps} {md_file}",
                       shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        result = subprocess.run(f"bash {scorecom} {md_file}", shell=True, check=True,
                                stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        it_score = float(result.stdout.strip())
        return i, it_score, md_file

# create DataFrame for an assembly path, following the same format as interfaces
def create_path_df(ifs, path):
    return pd.concat([ifs.iloc[0:0], pd.DataFrame(path)], ignore_index=True)

# sort transformations in the specified order
def sort_dataframe(df, score1, score2, chain1, chain2):
    data = df.to_numpy()
    score1_idx = df.columns.get_loc(score1)
    score2_idx = df.columns.get_loc(score2)
    chain1_idx = df.columns.get_loc(chain1)
    chain2_idx = df.columns.get_loc(chain2)
    sorted_idx = np.lexsort((
        data[:, chain2_idx],
        data[:, chain1_idx],
        -data[:, score2_idx],
        -data[:, score1_idx]
    ))
    sorted_data = data[sorted_idx]
    sorted_df = pd.DataFrame(sorted_data, columns=df.columns)
    return sorted_df

# randomly sort transformations
def shuffle_dataframe(df, seed, frac=1.0):
    n_rows = df.shape[0]
    sampled_size = int(frac * n_rows)
    np.random.seed(seed)
    shuffled_indices = np.random.permutation(n_rows)
    selected_indices = shuffled_indices[:sampled_size]
    shuffled_data = df.to_numpy()[selected_indices]
    return pd.DataFrame(shuffled_data, columns=df.columns)

# check the connectivity between two subunits derived from a long chain
def check_subunits_connectivity(subunit_id1, subunit_id2, subunit1, subunit2):
    original_chain_id = config.subunits_id_mapping[subunit_id1]
    start_res1, end_res1 = config.subunits_range[original_chain_id][subunit_id1]
    start_res2, end_res2 = config.subunits_range[original_chain_id][subunit_id2]
    if start_res2 > end_res1:
        # dist_threshold = abs(start_res2 - end_res1) * 4.0
        dist_threshold = (abs(start_res2 - end_res1) + 1) * 4.0
        coords1 = subunit1.residues[-1]["CA"].get_coord()
        coords2 = subunit2.residues[0]["CA"].get_coord()
    else:
        # dist_threshold = abs(start_res1 - end_res2) * 4.0
        dist_threshold = (abs(start_res1 - end_res2) + 1) * 4.0
        coords1 = subunit1.residues[0]["CA"].get_coord()
        coords2 = subunit2.residues[-1]["CA"].get_coord()
    breaking_distance = np.linalg.norm(coords1 - coords2)
    return breaking_distance <= dist_threshold

# convert split subunits into a full chain
def create_full_orientation(old_orientation):
    full_orientation = {}
    for chain_id in config.full_chain_ids:
        if chain_id not in old_orientation:
            continue
        if chain_id not in config.subunits_range:
            full_orientation[chain_id] = old_orientation[chain_id]
        else:
            full_chain = Chain(chain_id)
            subunit_ids = list(config.subunits_range[chain_id].keys())
            for subunit_id in subunit_ids:
                for residue in old_orientation[subunit_id].get_residues():
                    full_chain.residues.append(residue)
            full_orientation[chain_id] = full_chain
    return full_orientation

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
def calculate_rmsd_orientation(full_orientation1, full_orientation2):
    coords1 = np.array([
        residue["CA"].get_coord() for chain in full_orientation1.values()
        for residue in chain.get_residues() if "CA" in residue.atoms
    ])
    coords2 = np.array([
        residue["CA"].get_coord() for chain in full_orientation2.values()
        for residue in chain.get_residues() if "CA" in residue.atoms
    ])
    if coords1.shape != coords2.shape:
        raise ValueError("[ERROR] Not equal CA atoms number in compared models.")
    ca_rmsd = calculate_rmsd_coords(coords1, coords2)
    return ca_rmsd

# search for the optimal correspondence between subunits from two complexes
def search_subunits_mapping(full_orientation1, full_orientation2):
    centroids1 = {
        chain_id: calculate_centroid_chains(chain) for chain_id, chain in full_orientation1.items()
    }
    centroids2 = {
        chain_id: calculate_centroid_chains(chain) for chain_id, chain in full_orientation2.items()
    }
    best_mapping = {chain_id: chain_id for chain_id in config.full_chain_ids}
    max_iter = 10
    for _ in range(max_iter):
        improved = False
        for homo_chains in config.homo_chains_list:
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

    updated_full_orientation2 = {
        chain_id: full_orientation2[best_mapping[chain_id]] for chain_id in config.full_chain_ids
    }
    for chain_id in updated_full_orientation2:
        updated_full_orientation2[chain_id].id = chain_id
    return updated_full_orientation2

# calculate CA rmsd for two complexes
def calculate_rmsd(full_orientation1, full_orientation2):
    if config.homo_chains_list:
        full_orientation2 = search_subunits_mapping(full_orientation1, full_orientation2)
        ca_rmsd = calculate_rmsd_orientation(full_orientation1, full_orientation2)
        return ca_rmsd
    else:
        ca_rmsd = calculate_rmsd_orientation(full_orientation1, full_orientation2)
        return ca_rmsd

# calculate quality score for the population
def check_population_confidence(all_results, check_point, rmsd_threshold):
    all_results.sort(key=lambda x: x["if_score"], reverse=True)
    clstr_results = []
    for result1 in all_results:
        if "full_orientation" not in result1:
            full_orientation = create_full_orientation(result1["orientation"])
            result1["full_orientation"] = full_orientation
        added_to_clstr = False
        for result2 in clstr_results:
            ca_rmsd = calculate_rmsd(result1["full_orientation"], result2["full_orientation"])
            if ca_rmsd < rmsd_threshold:
                added_to_clstr = True
                break
        if not added_to_clstr:
            clstr_results.append(result1)
        if len(clstr_results) >= check_point:
            break
    confidence_score = clstr_results[-1]["if_score"]
    return confidence_score

# rank structures according to the specified score then cluster
def cluster_results(initial_results, population_size, rmsd_threshold, rank_type):
    if rank_type == 1:
        initial_results.sort(key=lambda x: x["if_score"], reverse=True)
    elif rank_type == 2:
        initial_results.sort(key=lambda x: x["it_score"])
    else:
        initial_results.sort(key=lambda x: x["snew_score"], reverse=True)
    clstr_results = []
    for result1 in initial_results:
        if "full_orientation" not in result1:
            full_orientation = create_full_orientation(result1["orientation"])
            result1["full_orientation"] = full_orientation
        added_to_clstr = False
        for result2 in clstr_results:
            ca_rmsd = calculate_rmsd(result1["full_orientation"], result2["full_orientation"])
            if ca_rmsd < rmsd_threshold:
                added_to_clstr = True
                break
        if not added_to_clstr:
            clstr_results.append(result1)
        if len(clstr_results) >= population_size:
            break
    if len(clstr_results) > population_size:
        clstr_results = clstr_results[:population_size]
    return clstr_results

# calculate IT-score of the structure and save it in the same dictionary
def add_results_itscore(results, num_cpus, md_steps, fout):
    stime = time.time()
    results_to_itscore = [
        {"idx": i, "result": result} for i, result in enumerate(results) if "it_score" not in result
    ]
    if not results_to_itscore:
        return results

    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        if md_steps > 0:
            futures = [
                executor.submit(
                    calculate_itscore_md, item["result"]["full_orientation"], item["idx"], md_steps
                ) for item in results_to_itscore
            ]
        else:
            futures = [
                executor.submit(
                    calculate_itscore, item["result"]["full_orientation"], item["idx"]
                ) for item in results_to_itscore
            ]
        for future in as_completed(futures):
            try:
                i, it_score, file = future.result()
                results[i]["it_score"] = it_score
                results[i]["file"] = file
            except Exception as e:
                print(f"[ERROR] Failed in calculating IT-score: {e}", flush=True)
            finally:
                del future
    gc.collect()
    etime = time.time()
    print(f"\ncalulate IT-score cost {etime - stime}", file=fout)
    return results

# refine the final output model, not used
def refine_output_model(i, result, model_file, md_steps):
    if md_steps > 0:
        file1 = f"{current_dir}/0-assemble_{i + 1}.pdb"
        write_pdb_file(result["full_orientation"], file1, config.format_lines)
        subprocess.run(f"python {openmm_script} {file1} {md_steps} {model_file}",
                       shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        os.remove(file1)
    else:
        write_pdb_file(result["full_orientation"], model_file, config.format_lines)

# print the final output model
def print_results_itscore(results, output, model_num, iteration, fout):
    results.sort(key=lambda x: x["it_score"])
    stem, suffix = output.rsplit(".", 1)
    num_models = len(results)
    num_models = min(num_models, model_num)
    selected_results = results[:num_models]
    for i, result in enumerate(selected_results):
        it_score = result["it_score"]
        print(f"Iteration {iteration} model {i + 1} it_score: {it_score}", file=fout)
        modle_file = f"{current_dir}/{stem}{iteration}_{i + 1}.{suffix}"
        shutil.copy(result["file"], modle_file)

# select a set of optimal subcomplexes when assembly fails
def select_optimal_subcomplex(initial_results):
    return max(
        initial_results,
        key=lambda x: (x["max_length"], x["mean_length"], x["if_score"])
    )

# comvert the group ID of the optimal subcomplexes into full chain ID
def create_full_group(old_groups):
    new_groups = {}
    for key, group in old_groups.items():
        mapped = set(
            config.subunits_id_mapping.get(chain_id, chain_id) for chain_id in group
        )
        unique_full = list(mapped)
        for chain_id in unique_full:
            new_groups[chain_id] = unique_full
    return new_groups

# filter chain based on the specified residue range
def filter_chain(chain, start_idx, end_idx):
    filtered_chain = Chain(chain.id)
    for residue in chain.get_residues():
        if residue.resseq in range(start_idx, end_idx + 1):
            filtered_chain.residues.append(residue)
    return filtered_chain

# print optimal subcomplexes for docking or subsequent processing
def print_groups(result, connectivity, fout):
    subcomplex_folder = f"{current_dir}/dock/"
    Path(subcomplex_folder).mkdir(exist_ok=True, parents=True)

    if not config.subunits_range:
        filtered_res_indices = config.filtered_res_indices
        groups, orientation = result["groups"], result["orientation"]
        all_chain_ids = config.chain_ids
    else:
        if connectivity:
            filtered_res_indices = {}
            for chain_id in config.full_chain_ids:
                if chain_id in config.subunits_range:
                    start_res_indices, end_res_indices = [], []
                    for subunit_id in config.subunits_range[chain_id]:
                        start_idx, end_idx = config.filtered_res_indices[subunit_id]
                        start_res_indices.append(start_idx)
                        end_res_indices.append(end_idx)
                    start_idx = min(start_res_indices)
                    end_idx = max(end_res_indices)
                    filtered_res_indices[chain_id] = (start_idx, end_idx)
                else:
                    filtered_res_indices[chain_id] = config.filtered_res_indices[chain_id]
            groups = create_full_group(result["groups"])
            orientation = create_full_orientation(result["orientation"])
            all_chain_ids = config.full_chain_ids
        else:
            filtered_res_indices = config.filtered_res_indices
            groups, orientation = result["groups"], result["orientation"]
            all_chain_ids = config.chain_ids
            with open(f"{subcomplex_folder}/warns.log", 'w+') as fwarn:
                print("[WARN] Split subunits fail to connect, please check manually!", file=fwarn)
            print("[WARN] Split subunits fail to connect, please check manually!", flush=True)

    in_group = set()
    subcomplexes = []
    for key, group in groups.items():
        if key not in in_group:
            new_group = ""
            for chain_id in group:
                in_group.add(chain_id)
                new_group += chain_id
            name = "".join(sorted(new_group))
            subcomplexes.append({"name": name, "model_type": 1})
    for chain_id in all_chain_ids:
        if chain_id not in in_group:
            in_group.add(chain_id)
            subcomplexes.append({"name": chain_id, "model_type": 0})

    subcomplexes.sort(key=lambda x: x["name"])
    for subcomplex in subcomplexes:
        name, model_type = subcomplex["name"], subcomplex["model_type"]
        model_file = f"{subcomplex_folder}/{name}.pdb"
        if model_type > 0:
            filtered_chains = []
            for chain_id in name:
                chain = orientation[chain_id]
                start_idx, end_idx = filtered_res_indices[chain_id]
                filtered_chain = filter_chain(chain, start_idx, end_idx)
                filtered_chains.append(filtered_chain)
            write_pdb_file(filtered_chains, model_file, config.format_lines)
        else:
            chain = config.mono_structs[name]
            start_idx, end_idx = filtered_res_indices[name]
            filtered_chain = filter_chain(chain, start_idx, end_idx)
            write_pdb_file(filtered_chain, model_file, config.format_lines)
    print("\nAssembly failed.", flush=True)
    print("Generated subcomplexes:", subcomplexes, flush=True)
    print("Generated subcomplexes:", subcomplexes, file=fout, flush=True)
