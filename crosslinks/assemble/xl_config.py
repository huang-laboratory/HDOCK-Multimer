# extract and prepare data for the assembly strategy

import pandas as pd
import numpy as np
import json
import tempfile
import os
import subprocess
import gc
from datetime import datetime
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from itertools import combinations
import sys

CONFIG_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(CONFIG_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
SCRIPTS_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "scripts"))
MMalign = os.path.join(TOOL_PATH, "MMalign")

sys.path.append(SCRIPTS_PATH)
from parse_pdb import PDBParser, load_format_line, write_pdb_file, Chain


pdb_parser = PDBParser()
RANDOM_SEED = 114514
steric_clash_dist = 2
BACKBONE = ("CA", "C", "N", "O")


chain_ids = None
chain_num = None
full_chain_ids = None
full_chain_num = None
mbfactors = None
mono_structs = None
source_structs = None
transformations = None
filtered_res_indices = None
homo_chains_list = None
total_length = None
subunits_range = None
subunits_order = None
subunits_id_mapping = None

chains_range = None
xl_restraints = None
pair_restraint_ids = None
all_xl_ids = None
xl_total_data = None

current_dir = os.getcwd()
temp_dir = os.path.join(current_dir, "temp")
os.makedirs(temp_dir, exist_ok=True)
format_lines = load_format_line(os.path.join(SCRIPTS_PATH, "format.pdb"))
interfaces = pd.read_csv(f"{current_dir}/interfaces.csv")


def get_num_cpus(max_workers: int):
    cpu_cap = int(os.getenv("NUM_CPUS", max_workers))
    cpu_count = os.cpu_count()
    return max(1, min(cpu_count, cpu_cap))

def filter_bfactor_all(chain: Chain, median_bfactor, min_bfactor):
    bfactor_threshold = min(median_bfactor, min_bfactor)
    filtered_chain = Chain(chain.id)
    residues = (
        residue for residue in chain.get_residues()
        if "CA" in residue.atoms and residue["CA"].get_bfactor() > bfactor_threshold
    )
    for residue in residues:
        filtered_chain.add_residue(residue)
    return filtered_chain

@contextmanager
def managed_tempfile(suffix="", delete_ok: bool = True):
    with tempfile.NamedTemporaryFile(
            suffix=suffix, dir=temp_dir, delete=delete_ok
    ) as temp_file:
        try:
            yield temp_file.name
        finally:
            temp_file.close()

def align(chain_id, source, fixed, moving, median_bfactor, min_bfactor):
    with managed_tempfile(".pdb") as fixed_file, \
            managed_tempfile(".pdb") as moving_file, \
            managed_tempfile() as mat_file:
        write_pdb_file(fixed, fixed_file, format_lines)
        if median_bfactor > 90:
            filtered_moving = filter_bfactor_all(moving, median_bfactor, min_bfactor)
            write_pdb_file(filtered_moving, moving_file, format_lines)
        else:
            write_pdb_file(moving, moving_file, format_lines)

        subprocess.run(f"{MMalign} {moving_file} {fixed_file} -m {mat_file} > /dev/null",
                       shell=True, check=True)
        rot, trans = [], []
        with open(mat_file, 'r') as file:
            next(file)
            next(file)
            for _ in range(3):
                line = file.readline().split()
                trans.append(float(line[1]))
                rot.append([float(i) for i in line[2:]])
        return chain_id, source, np.array(rot).T, np.array(trans)

def extract_median_bfactor(chain: Chain):
    plddt_values = [res["CA"].get_bfactor() for res in chain.get_residues() if "CA" in res.atoms]
    return np.median(plddt_values) if plddt_values else 0.0


def load_struct_data(pdbdir, file_stoi, num_cpus, fout):
    global chain_ids, chain_num, mbfactors, mono_structs, source_structs
    global transformations, filtered_res_indices, homo_chains_list
    global full_chain_ids, full_chain_num, total_length
    global subunits_range, subunits_order, subunits_id_mapping
    global chains_range

    unique_chains = pd.unique(interfaces[["chain1", "chain2"]].values.ravel("K"))
    chain_ids = sorted(unique_chains.tolist())
    chain_num = len(chain_ids)
    mbfactors = {}
    mono_structs = {}
    filtered_res_indices = {}
    homo_chains_list = []
    total_length = 0

    for chain_id in chain_ids:
        chain_struct = pdb_parser.get_structure("", f"{current_dir}/{chain_id}.pdb")[0][chain_id]
        median_bfactor = extract_median_bfactor(chain_struct)
        mbfactors[chain_id] = median_bfactor
        mono_structs[chain_id] = chain_struct

    with open(f"{current_dir}/res_indices_A.txt", 'r') as fidx:
        for line in fidx:
            parts = line.strip().split()
            filtered_res_indices[parts[0]] = (int(parts[1]), int(parts[2]))

    subunits_range, subunits_order, subunits_id_mapping = {}, {}, {}
    chains_range = {}
    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)
    unique_chain_ids = []
    for i in stoi_data:
        asym_ids, length = i["asym_ids"], i["length"]
        total_length += len(asym_ids) * length
        unique_chain_ids.extend(asym_ids)
        if len(asym_ids) > 1:
            homo_chains_list.append(asym_ids)
        if "subunits" in i:
            for chain_idx, chain_id in enumerate(asym_ids):
                subunits_range[chain_id], subunits_order[chain_id] = {}, {}
                for j in i["subunits"]:
                    subunit_order = int(j["subunit_id"].split("-")[-1])
                    subunit_id = j["asym_ids"][chain_idx]
                    start_res, end_res = j["start_res"], j["end_res"]
                    subunits_range[chain_id][subunit_id] = (start_res, end_res)
                    chains_range[subunit_id] = (start_res, end_res)
                    subunits_order[chain_id][subunit_id] = subunit_order
                    subunits_id_mapping[subunit_id] = chain_id
        else:
            for chain_id in asym_ids:
                chains_range[chain_id] = (1, length)

    full_chain_ids = sorted(unique_chain_ids)
    full_chain_num = len(full_chain_ids)

    stime = datetime.now()
    source_structs = {
        source: pdb_parser.get_structure("", os.path.join(pdbdir, source))[0]
        for source in interfaces["source"].unique()
    }
    transformations = {chain_id: {} for chain_id in chain_ids}
    futures = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        for chain_id in chain_ids:
            for source, source_struct in source_structs.items():
                if chain_id in source_struct.chains:
                    mono_chain = mono_structs[chain_id]
                    target_chain = source_struct[chain_id]
                    median_bfactor = mbfactors[chain_id]
                    futures.append(
                        executor.submit(
                            align, chain_id, source, target_chain, mono_chain, median_bfactor, 80
                        )
                    )
        for future in as_completed(futures):
            chain_id, source, rot, trans = future.result()
            transformations[chain_id][source] = (rot, trans)
    source_structs = {}
    del futures
    gc.collect()
    etime = datetime.now()
    print(f"Extracting transformations cost {etime - stime}\n", file=fout, flush=True)
    print(f"\nExtracting transformations cost {etime - stime}", flush=True)

def generate_chain_pair(chain1, chain2):
    return (chain1, chain2) if chain1 <= chain2 else (chain2, chain1)

def find_closest_res(res_seq, start_res, end_res, max_offset: int = 20):
    if start_res <= res_seq <= end_res:
        return res_seq, 0
    elif res_seq < start_res and start_res - res_seq <= max_offset:
        return start_res, start_res - res_seq
    elif res_seq > end_res and res_seq - end_res <= max_offset:
        return end_res, res_seq - end_res
    else:
        return None, None

def generate_xl_subunits(chains, res_seq, max_offset: int = 20):
    global subunits_range, subunits_id_mapping
    global chains_range

    xl_subunits = []
    for chain in set(chains):
        if subunits_range:
            if chain in subunits_range:
                candidates, exact_hits = [], []
                for subunit, (start_res, end_res) in subunits_range[chain].items():
                    hit = find_closest_res(res_seq, start_res, end_res, max_offset)
                    if not hit:
                        continue
                    xl_res, xl_offset = hit
                    entry = {"chain": subunit, "res_seq": xl_res, "offset": xl_offset}
                    candidates.append(entry)
                    if xl_offset == 0:
                        exact_hits.append(entry)
                if exact_hits:
                    xl_subunits.extend(exact_hits)
                else:
                    xl_subunits.extend(candidates)
            else:
                if chain not in chains_range:
                    continue
                start_res, end_res = chains_range[chain]
                hit = find_closest_res(res_seq, start_res, end_res, max_offset)
                if not hit:
                    continue
                xl_res, xl_offset = hit
                entry = {"chain": chain, "res_seq": xl_res, "offset": xl_offset}
                xl_subunits.append(entry)
        else:
            if chain not in chains_range:
                continue
            start_res, end_res = chains_range[chain]
            hit = find_closest_res(res_seq, start_res, end_res, max_offset)
            if not hit:
                continue
            xl_res, xl_offset = hit
            entry = {"chain": chain, "res_seq": xl_res, "offset": xl_offset}
            xl_subunits.append(entry)
    return xl_subunits

def load_crosslink_data(file_crosslink, max_offset: int = 20, allow_self_linked: bool = False):
    global xl_restraints, pair_restraint_ids, all_xl_ids, xl_total_data
    global chains_range, mono_structs, chain_ids

    xl_restraints = {}
    pair_restraint_ids = defaultdict(list)
    with open(file_crosslink, 'r') as flink:
        for line in flink:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue

            res_seq1, chains1 = int(parts[0]), parts[1]
            res_seq2, chains2 = int(parts[2]), parts[3]
            min_distance, max_distance = float(parts[4]), float(parts[5])
            w1 = float(parts[-1])

            xl_subunits1 = generate_xl_subunits(chains1, res_seq1, max_offset)
            if not xl_subunits1:
                continue
            xl_subunits2 = generate_xl_subunits(chains2, res_seq2, max_offset)
            if not xl_subunits2:
                continue

            chain_pairs = set()
            xl_id = len(xl_restraints)
            xl_derived_by_pairs = {}
            for sub1 in xl_subunits1:
                for sub2 in xl_subunits2:
                    chain1, chain2 = sub1["chain"], sub2["chain"]
                    if not allow_self_linked and chain1 == chain2:
                        continue
                    if chain1 not in chain_ids or chain2 not in chain_ids:
                        print(f"[WARN] Cannot find corresponding subunits for ({chain1}, {chain2})")
                        continue
                    pair = generate_chain_pair(chain1, chain2)
                    sub_res1, sub_res2 = sub1["res_seq"], sub2["res_seq"]
                    offset1, offset2 = sub1["offset"], sub2["offset"]
                    maxD_eff = max_distance + (offset1 + offset2) * 4.0
                    start_res1, end_res1 = chains_range[chain1]
                    start_res2, end_res2 = chains_range[chain2]
                    plddt1 = mono_structs[chain1].residues[sub_res1 - start_res1]["CA"].get_bfactor()
                    plddt2 = mono_structs[chain2].residues[sub_res2 - start_res2]["CA"].get_bfactor()
                    w2 = (plddt1 + plddt2) / 2.0
                    xl_derived_by_pairs[pair] = {
                        "chain1": chain1, "chain2": chain2, "res1": sub_res1, "res2": sub_res2,
                        "minD": min_distance, "maxD": maxD_eff, "w1": w1, "w2": w2
                    }
                    chain_pairs.add(pair)
            if not chain_pairs:
                continue
            xl_restraints[xl_id] = {
                "chains1": list(chains1), "chains2": list(chains2),
                "res1": res_seq1, "res2": res_seq2,
                "minD": min_distance, "maxD": max_distance, "w": w1,
                "pairs": sorted(chain_pairs),
                "xl_pairs": xl_derived_by_pairs
            }
            for pair in chain_pairs:
                pair_restraint_ids[pair].append(xl_id)

    all_pairs = [generate_chain_pair(c1, c2) for c1, c2 in combinations(chain_ids, 2)]
    all_xl_ids = set()
    for pair in all_pairs:
        if pair in pair_restraint_ids:
            all_xl_ids.update(set(pair_restraint_ids[pair]))
    xl_total_data = {}
    for xl_id in all_xl_ids:
        xl_pairs_data = xl_restraints[xl_id]["xl_pairs"]
        pair, pair_data = max(xl_pairs_data.items(), key=lambda item: item[1]["w2"])
        xl_total_data[xl_id] = {
            "pair": pair, "w1": pair_data["w1"], "w2": pair_data["w2"]
        }
