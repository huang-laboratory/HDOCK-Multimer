import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
import os
from scipy.spatial.distance import cdist
from itertools import combinations
from argparse import ArgumentParser
import json
import math
import time
import tempfile
from contextlib import contextmanager
import sys

DATA_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(DATA_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
SCRIPTS_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "scripts"))
STRIDE = os.path.join(TOOL_PATH, "stride")

sys.path.append(SCRIPTS_PATH)
from parse_pdb import PDBParser, load_format_line, write_pdb_file


protein_letters = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}
interaction_dist = 8
pdb_parser = PDBParser()
format_lines = load_format_line(os.path.join(SCRIPTS_PATH, "format.pdb"))
secondary_struct_flags = ("H", "G", "I", "B", "E")


# basic utility functions
@contextmanager
def managed_tempfile(suffix=""):
    with tempfile.NamedTemporaryFile(suffix=suffix, dir=os.getcwd()) as temp_file:
        try:
            yield temp_file.name
        finally:
            del temp_file

def calculate_mean_plddt(chain):
    plddt_values = [
        residue["CA"].get_bfactor() for residue in chain.get_residues() if "CA" in residue.atoms
    ]
    return np.mean(plddt_values)

def extract_secondary_structure(pdb_file, chain_id):
    secondary_structure = {}
    with managed_tempfile() as tempout:
        subprocess.run(f"{STRIDE} {pdb_file} > {tempout}", shell=True, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        with open(tempout, 'r') as file:
            for line in file:
                if line.startswith("ASG"):
                    parts = line.strip().split()
                    c_id, res_id, ss = parts[2], int(parts[3]), parts[5]
                    secondary_structure[(c_id, res_id)] = ss
    chain = pdb_parser.get_structure("", pdb_file)[0][chain_id]
    for residue in chain.get_residues():
        if (chain_id, residue.resseq) not in secondary_structure:
            secondary_structure[(chain_id, residue.resseq)] = "-"
    return secondary_structure

def extract_domains(chain, plddt_threshold, sstructure):
    domains = []
    current_domain = []
    for residue in chain.get_residues():
        res_id = residue.resseq
        ca_plddt = residue["CA"].get_bfactor()
        if ca_plddt >= plddt_threshold and sstructure[(chain.id, res_id)] in secondary_struct_flags:
            current_domain.append(res_id)
        else:
            if current_domain:
                domains.append(current_domain)
                current_domain = []
    if current_domain:
        domains.append(current_domain)

    filtered_domains = []
    for i, domain in enumerate(domains):
        if i > 0:
            dl1 = abs(domains[i - 1][-1] - domain[0])
        else:
            dl1 = float("inf")
        if i < len(domains) - 1:
            dl2 = abs(domain[-1] - domains[i + 1][0])
        else:
            dl2 = float("inf")
        dl = min(dl1, dl2)
        if len(domain) >= 10:
            filtered_domains.append(domain)
        elif len(domain) >= 4 and dl < 10:
            filtered_domains.append(domain)
    return filtered_domains

def filter_residue_ids(chain, mean_plddt, sstructure, keep_n=False, keep_c=False):
    if mean_plddt >= 80:
        plddt_threshold = 70
    elif mean_plddt >= 50:
        plddt_threshold = 50
    else:
        plddt_threshold = 0

    residues = list(chain.get_residues())
    if plddt_threshold >= 70:
        res_indices_N, res_indices_C = set(), set()
        keep_residues_N, keep_residues_C = keep_n, keep_c
        for residue in residues:
            res_id = residue.resseq
            ca_plddt = residue["CA"].get_bfactor()
            if keep_residues_N or ca_plddt >= plddt_threshold:
                res_indices_N.add(res_id)
                keep_residues_N = True
        for residue in reversed(residues):
            res_id = residue.resseq
            ca_plddt = residue["CA"].get_bfactor()
            if keep_residues_C or ca_plddt >= plddt_threshold:
                res_indices_C.add(res_id)
                keep_residues_C = True
        filtered_residue_indices = res_indices_N & res_indices_C
    elif plddt_threshold == 50:
        domains = extract_domains(chain, plddt_threshold, sstructure)
        if domains:
            min_idx = domains[0][0]
            max_idx = domains[-1][-1]
        else:
            min_idx = residues[0].resseq
            max_idx = residues[-1].resseq
        if keep_n:
            min_idx = residues[0].resseq
        if keep_c:
            max_idx = residues[-1].resseq
        filtered_residue_indices = set(range(min_idx, max_idx + 1))
    else:
        filtered_residue_indices = set(r.resseq for r in residues)
    return sorted(filtered_residue_indices)

def renumber_residue_ids(chain, start_res):
    for i, residue in enumerate(chain.residues):
        residue.resseq = i + start_res
    return chain


# Data Processor for assembly
class DataProcessor:
    def __init__(self, fout, pdbdir: Path, file_stoi: Path):
        self.fout = fout
        self.pdbdir = pdbdir
        self.stoi = file_stoi
        self.cwd = Path.cwd()

        self.chain_ids = []
        self.chain_num = 0
        self.sequences = {}
        self.homo_chains = {}
        self.mono_files = {}
        self.mean_plddt_values = {}
        self.start_res_indices = {}
        self.keep_chain_terminal = {}
        self.filtered_res_ids = {}

        self.read_stoi()

    # read stoichiometry json file
    def read_stoi(self):
        fastas = []
        with open(self.stoi, 'r') as fstoi:
            stoi_data = json.load(fstoi)

        for i in stoi_data:
            if "subunits" in i:
                length = i["length"]
                for j in i["subunits"]:
                    split_asym_ids, sequence = j["asym_ids"], j["sequence"]
                    start_res, end_res = j["start_res"], j["end_res"]
                    fastas.append(sequence)
                    seq = fastas.index(sequence)
                    self.chain_ids.extend(split_asym_ids)
                    for chain_id in split_asym_ids:
                        self.sequences[chain_id] = seq
                        self.start_res_indices[chain_id] = start_res
                        if start_res == 1:
                            self.keep_chain_terminal[chain_id] = (False, True)
                        elif end_res == length:
                            self.keep_chain_terminal[chain_id] = (True, False)
                        else:
                            self.keep_chain_terminal[chain_id] = (True, True)
                        if len(split_asym_ids) >= 5:
                            self.homo_chains[chain_id] = split_asym_ids
            else:
                asym_ids, sequence = i["asym_ids"], i["sequence"]
                fastas.append(sequence)
                seq = fastas.index(sequence)
                self.chain_ids.extend(asym_ids)
                for chain_id in asym_ids:
                    self.sequences[chain_id] = seq
                    self.start_res_indices[chain_id] = 1
                    self.keep_chain_terminal[chain_id] = (False, False)
                    if len(asym_ids) >= 5:
                        self.homo_chains[chain_id] = asym_ids
        self.chain_ids = sorted(self.chain_ids)
        self.chain_num = len(self.chain_ids)

    # judge whether the subcomponent is homologous
    def check_homology(self, chains_ids):
        homology = True
        chain_id1 = chains_ids[0]
        if chain_id1 not in self.homo_chains:
            homology = False
        else:
            for chain_id in chains_ids:
                if chain_id not in self.homo_chains[chain_id1]:
                    homology = False
                    break
        return homology

    # extract monomer and interface data from predicted subcomponents
    def extract_data_from_subcomplex(self):
        interfaces = {
            "source": [], "chain1": [], "chain2": [], "length1": [], "length2": [],
            "contacts": [], "score": [], "mean_score": [],
            "scaled_score1": [], "scaled_score2": []
        }
        monomers = {"source": [], "chain": [], "mean_plddt": []}
        seen_interfaces = {}
        stime = time.time()

        for subcomplex in self.pdbdir.glob("*.pdb"):
            structure = pdb_parser.get_structure("subcomplex", subcomplex)
            chains = list(structure.get_chains())
            source = subcomplex.name
            chains_residues = {chain.id: list(chain.get_residues()) for chain in chains}
            chains_ids = [chain.id for chain in chains]
            chains_seqs = {chain_id: self.sequences[chain_id] for chain_id in chains_ids}
            homology = self.check_homology(chains_ids)
            for chain1, chain2 in combinations(chains, 2):
                chain1_res = chains_residues[chain1.id]
                chain2_res = chains_residues[chain2.id]
                chain1_ca = np.array([res["CA"].get_coord() for res in chain1_res if "CA" in res.atoms])
                chain2_ca = np.array([res["CA"].get_coord() for res in chain2_res if "CA" in res.atoms])
                contacts = np.argwhere(cdist(chain1_ca, chain2_ca) < interaction_dist)
                if len(contacts) == 0:
                    continue

                interface_key = (chain1.id, chain2.id)
                other_seqs = tuple(
                    sorted(
                        (chains_seqs[chain.id] for chain in chains if chain.id not in interface_key)
                    )
                )
                chain1_inter, chain2_inter = set(contacts[:, 0]), set(contacts[:, 1])
                b_factors = np.array([chain1_res[i]["CA"].get_bfactor() for i in chain1_inter] +
                                     [chain2_res[i]["CA"].get_bfactor() for i in chain2_inter])
                inter_score = np.sum(b_factors)
                interface_avail = False
                if interface_key in seen_interfaces:
                    used_ifs = seen_interfaces[interface_key]
                    source_used = False
                    for used_source, used_other_seqs, used_score in used_ifs:
                        if used_other_seqs == other_seqs and abs(used_score - inter_score) <= 1e-6:
                            source_used = True
                            break
                    if not source_used:
                        seen_interfaces[interface_key].append((source, other_seqs, inter_score))
                        interface_avail = True
                else:
                    seen_interfaces[interface_key] = [(source, other_seqs, inter_score)]
                    interface_avail = True

                if interface_avail:
                    chain1_plddt = np.mean([res["CA"].get_bfactor() for res in chain1_res])
                    chain2_plddt = np.mean([res["CA"].get_bfactor() for res in chain2_res])
                    monomers["source"].append(source)
                    monomers["chain"].append(chain1.id)
                    monomers["mean_plddt"].append(chain1_plddt)
                    monomers["source"].append(source)
                    monomers["chain"].append(chain2.id)
                    monomers["mean_plddt"].append(chain2_plddt)

                    interfaces["source"].append(source)
                    chain1_id = min(chain1.id, chain2.id)
                    chain2_id = max(chain1.id, chain2.id)
                    interfaces["chain1"].append(chain1_id)
                    interfaces["chain2"].append(chain2_id)
                    interfaces["length1"].append(len(chains_residues[chain1_id]))
                    interfaces["length2"].append(len(chains_residues[chain2_id]))
                    contacts = len(chain1_inter) + len(chain2_inter)
                    mean_score = inter_score / contacts
                    interfaces["contacts"].append(contacts)
                    interfaces["score"].append(inter_score)
                    interfaces["mean_score"].append(mean_score)
                    # scaled_score1 = contacts ** 0.44 * mean_score
                    scaled_score1 = math.log(contacts, 40) * mean_score
                    if homology:
                        scaled_score2 = mean_score + mean_score * (1 - mean_score / 100)
                    else:
                        scaled_score2 = mean_score
                    interfaces["scaled_score1"].append(scaled_score1)
                    interfaces["scaled_score2"].append(scaled_score2)

        best_monomer = (
            pd.DataFrame.from_records(monomers).sort_values(by="mean_plddt", ascending=False).groupby("chain").first()
        )
        chain_ids = best_monomer.index.tolist()
        chain_num = len(chain_ids)
        if chain_num < self.chain_num:
            print("Assembly stage terminated.")
            raise SystemExit(f"[ERROR] Insufficient chains: {chain_num}, native chains: {self.chain_num}")
        for index, row in best_monomer.iterrows():
            chain_id = index
            source, plddt = row["source"], row["mean_plddt"]
            print(f"{chain_id} from subcomplex {source}, mean plddt {plddt}", file=self.fout)
            chain_file = self.cwd / f"{chain_id}.pdb"
            self.mean_plddt_values[chain_id] = plddt
            self.mono_files[chain_id] = chain_file
            chain = pdb_parser.get_structure("", self.pdbdir / source)[0][chain_id]
            start_res = self.start_res_indices[chain_id]
            renumber_residue_ids(chain, start_res)
            write_pdb_file(chain, chain_file, format_lines)

        interfaces = pd.DataFrame.from_dict(interfaces)
        interfaces.sort_values(
            by=["chain1", "chain2", "mean_score"], ascending=[True, True, False], inplace=True
        )
        interfaces.to_csv(self.cwd / "interfaces.csv", index=False)
        etime = time.time()
        print(f"\nExtract monomers and interfaces cost {etime - stime}", file=self.fout, flush=True)
        print(f"\nExtract monomers and interfaces cost {etime - stime}", flush=True)

    # filter low-quality terminal of extracted representative monomers
    def filter_mono_chains(self):
        for chain_id in self.chain_ids:
            chain_file = self.mono_files[chain_id]
            chain = pdb_parser.get_structure("", chain_file)[0][chain_id]
            mean_plddt = calculate_mean_plddt(chain)
            sstructure = extract_secondary_structure(chain_file, chain_id)
            keep_n, keep_c = self.keep_chain_terminal[chain_id]
            filtered_residue_indices = filter_residue_ids(
                chain, mean_plddt, sstructure, keep_n, keep_c
            )
            self.filtered_res_ids[chain_id] = (
                filtered_residue_indices[0], filtered_residue_indices[-1]
            )


####################      MAIN      ####################
def main():
    arg_parser = ArgumentParser()
    arg_parser.add_argument("-d", "--pdbdir", type=Path, required=True)
    arg_parser.add_argument("-s", "--stoi", type=Path, required=True)

    args = arg_parser.parse_args()
    subcomplex_dir = args.pdbdir
    file_stoi = args.stoi

    current_dir = Path.cwd()
    with open(current_dir / "data.log", 'w+') as fout:
        data_processor = DataProcessor(fout, subcomplex_dir, file_stoi)
        data_processor.extract_data_from_subcomplex()
        data_processor.filter_mono_chains()

    with open(current_dir / "res_indices_A.txt", 'w+') as fidx:
        for chain_id in data_processor.chain_ids:
            start_idx, end_idx = data_processor.filtered_res_ids[chain_id]
            print(f"{chain_id}    {start_idx}    {end_idx}", file=fidx)


if __name__ == "__main__":
    main()

