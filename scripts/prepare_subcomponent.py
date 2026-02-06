<<<<<<< HEAD
# scripts for generating subcomponent structures for repeated subunits

=======
>>>>>>> af68c8e3e373774427ca682a31edb6f375c5b95e
import sys
import json
from pathlib import Path
from itertools import combinations, permutations, combinations_with_replacement, product
import os
import shutil


def generate_nmer(file_stoi, pdb_id: str, n: int, temp_dir: str):
    with open(file_stoi, 'r') as fstoi:
        stoi_data = json.load(fstoi)

    unique_chain_ids = {}
    entity_ids = {}
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
    chain_ids = sorted(list(unique_chain_ids.keys()))
    print(unique_chain_ids)

    pdb_id = pdb_id.upper()
    current_dir = os.path.dirname(os.path.abspath(file_stoi))
    subcomponent_dir = Path(f"{current_dir}/subcomponents/")
    subcomponent_dir.mkdir(exist_ok=True, parents=True)

    stoi = {}
    subcomponent_indice_file = f"{current_dir}/new_subcomponent.txt"
    with open(subcomponent_indice_file, 'w+') as fidx:
        for name in combinations_with_replacement(chain_ids, n):
            for chain_id in chain_ids:
                stoi[chain_id] = name.count(chain_id)

            subunit_indice = "_".join(entity_ids[chain_id] for chain_id in name)
            temp_file = Path(temp_dir) / Path(f"{pdb_id}_{subunit_indice}.pdb")
            if not temp_file.exists():
                continue

            struct = {chr(ord("A") + i): [] for i in range(n)}
            with open(temp_file, 'r') as file:
                for line in file:
                    if line.startswith("ATOM"):
                        struct[line[21]].append(line)

            permutation_num = [0 for _ in range(len(name))]
            for i, chain_id in enumerate(name):
                permutation_num[i] = stoi[chain_id]
                stoi[chain_id] = 0
            perms_all = []
            for idx in range(len(name)):
                chain_id = name[idx]
                count = permutation_num[idx]
                perms = list(permutations(unique_chain_ids[chain_id], count))
                perms_all.append(perms)
            for comb in product(*perms_all):
                chains = []
                for group in comb:
                    chains.extend(group)
                subcomponent_filename = pdb_id + "_" + "".join(chains) + ".pdb"
                subcomponent_file = subcomponent_dir / Path(subcomponent_filename)
                with open(subcomponent_file, 'w') as file:
                    for i in range(n):
                        chain_id = chains[i]
                        for line in struct[chr(ord("A") + i)]:
                            file.write(line[:21] + chain_id + line[22:])
                print(subcomponent_filename, file=fidx)


def main():
    file_stoi = sys.argv[1]
    pdb_id = sys.argv[2]
    n = int(sys.argv[3])
    temp_dir = sys.argv[4]

    generate_nmer(file_stoi, pdb_id, n, temp_dir)
    

if __name__ == "__main__":
    main()
