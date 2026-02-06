# modify the chain IDs in the docking models

import json
import sys
import os
from pathlib import Path

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
sys.path.append(BASE_PATH)

from parse_pdb import PDBParser, load_format_line, write_pdb_file
from parse_pdb import Chain

protein_letters = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}
pdb_parser = PDBParser()
format_lines = load_format_line(os.path.join(BASE_PATH, "format.pdb"))


# extract sequence for the given chain, "X" for nonstandard residues
def extract_sequence(chain: Chain):
    return "".join(protein_letters.get(res.resname.strip(), "X") for res in chain.get_residues())


############    MAIN     ############
def main():
    file_stoi = sys.argv[1]

    stoi_dir = os.path.dirname(os.path.abspath(file_stoi))
    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    file1 = f"{stoi_dir}/asym/models/multi.pdb"
    file2 = f"{stoi_dir}/asym/models/new_multi.pdb"

    model = pdb_parser.get_structure("", file1)[0]
    chains = list(model.get_chains())
    orientation = {chain.id: chain for chain in chains}
    chain_seqs = {chain.id: extract_sequence(chain) for chain in chains}
    unique_chain_ids = {}
    unique_chain_seqs = {}
    for i in stoi_data:
        asym_ids, sequence = i["asym_ids"], i["sequence"]
        unique_chain_ids[asym_ids[0]] = asym_ids
        unique_chain_seqs[asym_ids[0]] = sequence
    chain_ids_maps = {}
    for chain_id1 in orientation.keys():
        chain_seq = chain_seqs[chain_id1]
        seq_found = False
        for chain_id2 in unique_chain_ids.keys():
            unique_chain_seq = unique_chain_seqs[chain_id2]
            if chain_seq == unique_chain_seq or chain_seq in unique_chain_seq:
                if chain_id2 in chain_ids_maps:
                    chain_ids_maps[chain_id2].append(chain_id1)
                else:
                    chain_ids_maps[chain_id2] = [chain_id1]
                seq_found = True
                break
        if not seq_found:
            print(f"{chain_id1}:", chain_seq)
            raise ValueError(f"[ERROR] Wrong chain sequence for {file1}")

    unique_chain_id_maps = {}
    for chain_id in unique_chain_ids.keys():
        old_chain_ids = chain_ids_maps[chain_id]
        new_chain_ids = unique_chain_ids[chain_id]
        for i, old_chain_id in enumerate(old_chain_ids):
            unique_chain_id_maps[old_chain_id] = new_chain_ids[i]
    new_orientation = {}
    for old_chain_id, new_chain_id in unique_chain_id_maps.items():
        chain = orientation[old_chain_id].copy()
        chain.id = new_chain_id
        new_orientation[new_chain_id] = chain

    write_pdb_file(new_orientation, file2, format_lines)
    os.remove(file1)
    os.rename(file2, file1)


if __name__ == "__main__":
    main()
