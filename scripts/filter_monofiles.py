import numpy as np
import subprocess
from parse_pdb import PDBParser, load_format_line, write_pdb_file
from parse_pdb import Chain
import json
import sys
import os
import time
import tempfile
from contextlib import contextmanager


BASE_PATH = os.path.dirname(os.path.abspath(__file__))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
STRIDE = os.path.join(TOOL_PATH, "stride")
secondary_struct_flags = ("H", "G", "I", "B", "E")
protein_letters = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}
pdb_parser = PDBParser()
format_lines = load_format_line(os.path.join(BASE_PATH, "format.pdb"))


@contextmanager
def managed_tempfile(suffix=""):
    with tempfile.NamedTemporaryFile(suffix=suffix, dir=os.getcwd()) as temp_file:
        try:
            yield temp_file.name
        finally:
            del temp_file

def extract_sequence(chain):
    return "".join(protein_letters.get(residue.resname.strip(), "X") for residue in chain.residues)

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

def filter_residue_ids(chain, mean_plddt, sstructure):
    if mean_plddt >= 80:
        plddt_threshold = 70
    elif mean_plddt >= 50:
        plddt_threshold = 50
    else:
        plddt_threshold = 0

    residues = list(chain.get_residues())
    if plddt_threshold >= 70:
        res_indices_N, res_indices_C = set(), set()
        keep_residues_N, keep_residues_C = False, False
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
        filtered_residue_indices = set(range(min_idx, max_idx + 1))
    else:
        filtered_residue_indices = set(r.resseq for r in residues)

    filtered_residue_indices = sorted(filtered_residue_indices)
    return filtered_residue_indices


def main():
    if len(sys.argv) < 3:
        print("Usage: python filter_monofiles.py <stoichiometry.json> <mono_dir>")
        sys.exit(1)

    file_stoi = sys.argv[1]
    mono_dir = sys.argv[2]

    start_time = time.time()
    current_dir = os.path.dirname(os.path.abspath(file_stoi))
    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    unique_chain_ids = {}
    sequences = {}
    for i in stoi_data:
        asym_ids, sequence = i["asym_ids"], i["sequence"]
        unique_chain_ids[asym_ids[0]] = asym_ids
        for chain_id in asym_ids:
            sequences[chain_id] = sequence

    filtered_res_ids = {}
    for chain_id in unique_chain_ids.keys():
        file1 = os.path.join(mono_dir, f"{chain_id}.pdb")
        file2 = f"{current_dir}/{chain_id}.pdb"
        #chain = pdb_parser.get_structure("", file1)[0][chain_id]
        model = pdb_parser.get_structure("", file1)[0]
        if chain_id in model.chains:
            chain = model[chain_id]
            mono_chain_id = chain_id
        elif "A" in model.chains:
            chain = model["A"]
            mono_chain_id = "A"
        else:
            raise ValueError(f"[ERROR] Chain {chain_id} not found")
        sequence = extract_sequence(chain)
        if sequence != sequences[chain_id]:
            raise ValueError(f"[ERROR] Wrong sequence for monomer: {chain_id}")

        mean_plddt = calculate_mean_plddt(chain)
        sstructure = extract_secondary_structure(file1, mono_chain_id)
        filtered_residue_indices = filter_residue_ids(chain, mean_plddt, sstructure)
        filtered_chain = Chain(chain_id)
        for residue in chain.get_residues():
            if residue.resseq in filtered_residue_indices:
                filtered_chain.add_residue(residue)
        write_pdb_file(filtered_chain, file2, format_lines)
        #os.remove(file1)
        #os.rename(file2, file1)
        filtered_res_ids[chain_id] = (
            filtered_residue_indices[0], filtered_residue_indices[-1]
        )

    with open(f"{current_dir}/res_indices_D.txt", 'w+') as fidx:
        for chain_id in unique_chain_ids.keys():
            start_idx, end_idx = filtered_res_ids[chain_id]
            print(f"{chain_id}    {start_idx}    {end_idx}", file=fidx)

    end_time = time.time()
    print(f"\nFilter monomer files cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    main()
