#!/data/xyao/anaconda3/bin/python

from Bio.PDB import PDBParser
import sys
from argparse import ArgumentParser
import numpy as np
import random
import json
from scipy.spatial.distance import cdist
import os


pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)

def load_chains_seqs(file_stoi):
    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    chain_to_seq = {}
    seq_to_chain = {}
    ident_chain_map = {}
    for i in stoi_data:
        asym_ids, sequence = i["asym_ids"], i["sequence"]
        for chain_id in asym_ids:
            chain_to_seq[chain_id] = sequence
            if sequence in seq_to_chain:
                seq_to_chain[sequence].append(chain_id)
            else:
                seq_to_chain[sequence] = [chain_id]
            ident_chain_map[chain_id] = asym_ids[0]
    return chain_to_seq, seq_to_chain, ident_chain_map


def load_res_plddt(chain_ids, multimer_dir):
    res_plddts = {}
    for chain_id in chain_ids:
        chain_file = os.path.join(multimer_dir, f"{chain_id}.pdb")
        chain_struct = pdb_parser.get_structure("", chain_file)
        chain_model = next(iter(chain_struct))
        for residue in chain_model.get_residues():
            if "CA" not in residue:
                continue
            res_plddts[(residue.parent.id, residue.id[1])] = residue["CA"].get_bfactor()
    return res_plddts


def simulate_crosslinks(file_stoi, file_native, multimer_dir):
    chain_to_seq, seq_to_chain, ident_chain_map = load_chains_seqs(file_stoi)
    chain_ids = sorted(chain_to_seq)
    res_plddts = load_res_plddt(chain_ids, multimer_dir)
    # np.random.seed(0)
    random.seed(0)
    native_struct = pdb_parser.get_structure("native", file_native)
    native_model = next(iter(native_struct))

    coords, identifiers = [], []
    for residue in native_model.get_residues():
        if residue.get_resname() == "LYS" and "CA" in residue:
            coords.append(residue["CA"].get_coord())
            identifiers.append((residue.parent.id, residue.id[1]))

    coord_distances = cdist(coords, coords)
    for i in range(len(coord_distances)):
        coord_distances[i][i:] = np.inf

    close_residues = np.argwhere(coord_distances < 30)
    far_residues = np.argwhere(coord_distances > 35)
    inter_close_residues = [i for i in close_residues if identifiers[i[0]][0] != identifiers[i[1]][0]]
    inter_far_residues = [i for i in far_residues if identifiers[i[0]][0] != identifiers[i[1]][0]]
    print(f"Total {len(coords)} LYS in the complex", flush=True)
    print(f"Total {len(close_residues)} close (distance < 30 A) LYS pairs in the complex", flush=True)
    print(f"Total {len(inter_close_residues)} inter close LYS pairs in the complex", flush=True)

    filtered_res = []
    np_coords = np.array(coords)
    for res1, res2 in inter_close_residues:
        coord1, coord2 = coords[res1], coords[res2]
        dir_vec = coord2 - coord1
        dir_size = np.linalg.norm(dir_vec)
        disturbed = False
        for i in range(3, int(dir_size) - 3, 2):
            checked_coord = coord1 + dir_vec * (i / dir_size)
            disturb_res = np.argwhere(np.linalg.norm(np_coords - checked_coord, axis=1) < 1)
            if len(disturb_res) > 0:
                disturbed = True
                break
        if not disturbed:
            filtered_res.append((identifiers[res1], identifiers[res2]))
    print(f"There are {len(filtered_res)} pairs of undisturbed residues remained")

    random.seed(114514)
    crosslinks_plddt = {}
    for (chain1, resseq1), (chain2, resseq2) in filtered_res:
        plddt1 = res_plddts[(chain1, resseq1)]
        plddt2 = res_plddts[(chain2, resseq2)]
        w1 = round((plddt1 + plddt2) / 200, 2)
        rep_link = tuple(sorted(
            ((ident_chain_map[chain1], resseq1), (ident_chain_map[chain2], resseq2))
        ))
        if rep_link in crosslinks_plddt:
            crosslinks_plddt[rep_link] = max(crosslinks_plddt[rep_link], w1)
        else:
            crosslinks_plddt[rep_link] = w1
    crosslinks = sorted(crosslinks_plddt.keys())

    k = int(len(crosslinks) * 0.1)
    selected_crosslinks = set(random.sample(crosslinks, k)) if k > 0 else set()
    print(k)
    if not selected_crosslinks:
        raise ValueError("No true crosslinks selected, please check complex")

    false_crosslinks_plddt = {}
    for res1, res2 in inter_far_residues:
        chain1, resseq1 = identifiers[res1]
        chain2, resseq2 = identifiers[res2]
        plddt1 = res_plddts[(chain1, resseq1)]
        plddt2 = res_plddts[(chain2, resseq2)]
        w1 = round((plddt1 + plddt2) / 200, 2)
        rep_link = tuple(sorted(
            ((ident_chain_map[chain1], resseq1), (ident_chain_map[chain2], resseq2))
        ))
        if rep_link in false_crosslinks_plddt:
            false_crosslinks_plddt[rep_link] = max(false_crosslinks_plddt[rep_link], w1)
        else:
            false_crosslinks_plddt[rep_link] = w1
    false_crosslinks = sorted(false_crosslinks_plddt.keys())
    kf = int(len(selected_crosslinks) * 0.05)
    print(kf)
    selected_false_crosslinks = set(random.sample(false_crosslinks, kf)) if kf > 0 else set()

    total_crosslinks = sorted(selected_crosslinks | selected_false_crosslinks)
    xl_lines = []
    for (chain1, resseq1), (chain2, resseq2) in total_crosslinks:
        chains1 = [k for k, v in ident_chain_map.items() if v == chain1]
        chains2 = [k for k, v in ident_chain_map.items() if v == chain2]
        w1 = crosslinks_plddt.get(((chain1, resseq1), (chain2, resseq2)))
        if w1 is None:
            w1 = false_crosslinks_plddt[((chain1, resseq1), (chain2, resseq2))]
        xl_line = f"{resseq1} {''.join(chains1)} {resseq2} {''.join(chains2)} 0 30 {w1}\n"
        xl_lines.append(xl_line)

    return xl_lines


def main():
    arg_parser = ArgumentParser()
    arg_parser.add_argument("-s", "--stoi", type=str, required=True)
    arg_parser.add_argument("-p", "--pdb", type=str, required=True)
    arg_parser.add_argument("-d", "--pdbdir", type=str, required=True)
    arg_parser.add_argument("-o", "--output", type=str, required=True)

    args = arg_parser.parse_args()
    file_stoi = args.stoi
    file_native = args.pdb
    multimer_dir = args.pdbdir
    file_output = args.output

    simulated_xl_lines = simulate_crosslinks(file_stoi, file_native, multimer_dir)
    with open(file_output, 'w+') as fout:
        fout.writelines(simulated_xl_lines)


if __name__ == "__main__":
    main()
