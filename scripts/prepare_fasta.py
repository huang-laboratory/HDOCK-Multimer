import sys
import json
from typing import List, Dict
from itertools import combinations, permutations, combinations_with_replacement
import os
from argparse import ArgumentParser


def generate_subcomponents_indices(stoi_data, subunits: int):
    unique_chain_ids: Dict[str, List[str]] = {}
    entity_ids: Dict[str, str] = {}
    for i in stoi_data:
        if "subunits" in i:
            for j in i["subunits"]:
                split_asym_ids = j["asym_ids"]
                unique_chain_ids[split_asym_ids[0]] = split_asym_ids
                entity_ids[split_asym_ids[0]] = j["subunit_id"]
        else:
            asym_ids = i["asym_ids"]
            unique_chain_ids[asym_ids[0]] = asym_ids
            entity_ids[asym_ids[0]] = i["entity_id"]
    chain_ids = list(unique_chain_ids.keys())
    chain_ids = sorted(chain_ids)

    subcomponents_indices = []
    for name in combinations_with_replacement(chain_ids, subunits):
        exist_ok = True
        for chain_id in chain_ids:
            if name.count(chain_id) > len(unique_chain_ids[chain_id]):
                exist_ok = False
                break
        if not exist_ok:
            continue
        subcomponent = []
        for chain_id in name:
            subcomponent.append(entity_ids[chain_id])
        subcomponents_indices.append(subcomponent)
    return subcomponents_indices


def generate_mono_fastas(fasta_dir, stoi_data, pdb_id):
    pdb_id = pdb_id.upper()
    for i in stoi_data:
        entity_id, sequence = i["entity_id"], i["sequence"]
        length = i["length"]
        file_fasta = f"{fasta_dir}/{pdb_id}_{entity_id}.fasta"
        lines_fasta = [
            f">{pdb_id} | Chain {entity_id} | {length}aa\n",
            f"{sequence}\n"
        ]
        with open(file_fasta, 'w+') as fout:
            fout.writelines(lines_fasta)


def generate_subcomponent_fastas(subcomp_dir, stoi_data, indices, pdb_id):
    pdb_id = pdb_id.upper()
    entity_id_sequences: Dict[str, str] = {}
    for i in stoi_data:
        if "subunits" in i:
            for j in i["subunits"]:
                subunit_id, sequence = j["subunit_id"], j["sequence"]
                entity_id_sequences[subunit_id] = sequence
        else:
            entity_id, sequence = i["entity_id"], i["sequence"]
            entity_id_sequences[entity_id] = sequence

    for i in indices:
        subunit_ids = "_".join(i)
        file_fasta = f"{subcomp_dir}/{pdb_id}_{subunit_ids}.fasta"
        lines_fasta = []
        for j in i:
            sequence = entity_id_sequences[j]
            length = len(sequence)
            lines_fasta.append(f">{pdb_id} | Chain {j} | {length}aa\n")
            lines_fasta.append(f"{sequence}\n")
        with open(file_fasta, 'w+') as fout:
            fout.writelines(lines_fasta)


def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--stoi", type=str, required=True)
    parser.add_argument("-n", "--nmer", type=int, default=3)
    parser.add_argument("-p", "--pdb", type=str, required=True)

    args = parser.parse_args()
    file_stoi = args.stoi
    nmer = args.nmer
    pdb_id = args.pdb

    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    # create directory for fasta files
    stoi_dir = os.path.dirname(os.path.abspath(file_stoi))
    mono_fasta_dir = os.path.join(stoi_dir, "mono_fastas")
    subcomponent_fasta_dir = os.path.join(stoi_dir, "subcomp_fastas")
    os.makedirs(mono_fasta_dir, exist_ok=True)
    os.makedirs(subcomponent_fasta_dir, exist_ok=True)

    # generate fasta files for monomers (full chains)
    generate_mono_fastas(mono_fasta_dir, stoi_data, pdb_id)

    # generate fasta files for subcomponents (n-mers)
    nmer_indices = generate_subcomponents_indices(stoi_data, nmer)
    generate_subcomponent_fastas(
        subcomponent_fasta_dir, stoi_data, nmer_indices, pdb_id
    )


if __name__ == "__main__":
    main()
