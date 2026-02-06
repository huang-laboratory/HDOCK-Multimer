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
import sys

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
MMalign = os.path.join(TOOL_PATH, "MMalign")

sys.path.append(BASE_PATH)
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
subunits_range = None
subunits_order = None
subunits_id_mapping = None

current_dir = os.getcwd()
temp_dir = os.path.join(current_dir, "temp")
os.makedirs(temp_dir, exist_ok=True)
format_lines = load_format_line(os.path.join(BASE_PATH, "format.pdb"))
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
    global full_chain_ids, full_chain_num
    global subunits_range, subunits_order, subunits_id_mapping

    unique_chains = pd.unique(interfaces[["chain1", "chain2"]].values.ravel("K"))
    chain_ids = sorted(unique_chains.tolist())
    chain_num = len(chain_ids)
    mbfactors = {}
    mono_structs = {}
    filtered_res_indices = {}
    homo_chains_list = []

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
    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)
    unique_chain_ids = []
    for i in stoi_data:
        asym_ids = i["asym_ids"]
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
                    subunits_order[chain_id][subunit_id] = subunit_order
                    subunits_id_mapping[subunit_id] = chain_id
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
