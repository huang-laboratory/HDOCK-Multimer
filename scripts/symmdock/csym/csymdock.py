import os
import shutil
import subprocess
import gc
import sys
import glob
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, "..", ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))

sys.path.append(BASE_PATH)
from parse_pdb import PDBParser

scorecom = os.path.join(TOOL_PATH, "scorecom.sh")
chdock = os.path.join(TOOL_PATH, "chdock")
compcn = os.path.join(TOOL_PATH, "compcn")
splitmodels_script = os.path.join(TOOL_PATH, "splitmodels.py")
pdb_parser = PDBParser()


############    Utility Functions    ############
# select appropriate number of parallel CPU cores
def get_num_cpus(max_workers: int):
    cpu_cap = int(os.getenv("NUM_CPUS", max_workers))
    cpu_count = os.cpu_count()
    return max(1, min(cpu_count, cpu_cap))

# renumber the chains in the file starting from "A"
def change_chain_ids(input_file, output_file):
    chain_ids = ([chr(i) for i in range(ord("A"), ord("Z") + 1)] +
                 [chr(i) for i in range(ord("a"), ord("z") + 1)] +
                 [chr(i) for i in range(ord("0"), ord("9") + 1)])
    current_chain = None
    output_lines = []
    chain_idx = 0

    with open(input_file, 'r') as fin:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                if current_chain is None:
                    current_chain = chain_ids[chain_idx]
                new_line = line[:21] + current_chain + line[22:]
                output_lines.append(new_line)
            elif line.startswith("TER"):
                if current_chain is None:
                    current_chain = chain_ids[chain_idx]
                new_line = line[:21] + current_chain + line[22:]
                output_lines.append(new_line)

                chain_idx += 1
                if chain_idx > len(chain_ids):
                    print(f"[ERROR] Too many chains in input PDB or too few chain IDS provided.",
                          file=sys.stderr)
                    sys.exit(1)
                current_chain = chain_ids[chain_idx] if chain_idx < len(chain_ids) else None

    with open(output_file, 'w+') as fout:
        fout.writelines(output_lines)

# define the cyclic symmetric sampling process
def run_csym_dock(stem, idx, csym, nmax):
    stem_file = f"{stem}_{idx}.pdb"
    out_file = f"{stem}_{idx}_c{csym}.out"
    log_file = f"cn_{idx}.log"
    out_pdb_file = f"{stem}_{idx}_{csym}.pdb"

    if not Path(stem_file).exists():
        print(f"[WARN] {stem_file} not found, skipping.", file=sys.stderr)
        return 0
    subprocess.run(f"{chdock} {stem_file} {stem_file} -symmetry {csym} -out {out_file} > {log_file} &&"
                   f"{compcn} {out_file} {out_pdb_file} -nmax {nmax} -complex > /dev/null &&"
                   f"python {splitmodels_script} {out_pdb_file}", shell=True, check=True)
    os.remove(out_pdb_file)

    for j in range(1, nmax+1):
        file1 = f"{stem}_{idx}_{csym}_{j}.pdb"
        if not Path(file1).is_file():
            raise RuntimeError(f"[ERROR] Change chain IDs failed for {file1}.")
        model_idx = nmax * (idx - 1) + j
        file2 = f"model_{model_idx}.pdb"
        # subprocess.run(f"changecid.py {file1} {file2}")
        change_chain_ids(file1, file2)
        os.remove(file1)

# calculate IT-score for the given model
def calculate_itscore_model(model):
    model_file = model["file"]
    model_idx = model["idx"]
    result = subprocess.run(f"bash {scorecom} {model_file}", shell=True, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    it_score = float(result.stdout.strip())
    return model_idx, it_score

# calculate IT-score for all models
def calculate_itscores_all(models, num_cpus):
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [
            executor.submit(calculate_itscore_model, model)
            for model in models
        ]
        for future in as_completed(futures):
            i, it_score = future.result()
            models[i]["it_score"] = it_score
    gc.collect()
    return models

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

# calculate CA rmsd for two complexes
def calculate_rmsd_orientation(orientation1, orientation2):
    coords1 = np.array([
        residue["CA"].get_coord() for chain in orientation1.values()
        for residue in chain.get_residues() if "CA" in residue.atoms
    ])
    coords2 = np.array([
        residue["CA"].get_coord() for chain in orientation2.values()
        for residue in chain.get_residues() if "CA" in residue.atoms
    ])
    if coords1.shape != coords2.shape:
        raise ValueError("[ERROR] Not equal CA atoms number in compared models.")
    ca_rmsd = calculate_rmsd_coords(coords1, coords2)
    return ca_rmsd

# search for the optimal correspondence between subunits from two complexes
def search_subunits_mapping(orientation1, orientation2, homo_chains_list):
    centroids1 = {chain_id: calculate_centroid_chains(chain) for chain_id, chain in orientation1.items()}
    centroids2 = {chain_id: calculate_centroid_chains(chain) for chain_id, chain in orientation2.items()}
    best_mapping = {chain_id: chain_id for chain_id in homo_chains_list}  # Homomeric complex
    max_iter = 10
    for _ in range(max_iter):
        improved = False
        homo_best_mapping = {chain_id: best_mapping[chain_id] for chain_id in homo_chains_list}
        homo_centroids1 = np.array([centroids1[chain_id] for chain_id in homo_chains_list])
        homo_centroids2 = np.array([
            centroids2[homo_best_mapping[chain_id]] for chain_id in homo_chains_list
        ])
        best_rmsd = calculate_rmsd_coords(homo_centroids1, homo_centroids2)
        for i, chain1 in enumerate(homo_chains_list):
            for j, chain2 in enumerate(homo_chains_list):
                if i >= j:
                    continue
                swapped_mapping = homo_best_mapping.copy()
                swapped_mapping[chain1], swapped_mapping[chain2] = swapped_mapping[chain2], swapped_mapping[chain1]
                homo_centroids2 = np.array([
                    centroids2[swapped_mapping[chain_id]] for chain_id in homo_chains_list
                ])
                current_rmsd = calculate_rmsd_coords(homo_centroids1, homo_centroids2)
                if current_rmsd < best_rmsd:
                    best_rmsd = current_rmsd
                    homo_best_mapping = swapped_mapping
                    improved = True
        if improved:
            for chain_id in homo_chains_list:
                best_mapping[chain_id] = homo_best_mapping[chain_id]
        else:
            break

    updated_orientation2 = {
        chain_id: orientation2[best_mapping[chain_id]].copy() for chain_id in homo_chains_list
    }
    for chain_id in updated_orientation2:
        updated_orientation2[chain_id].id = chain_id
    return updated_orientation2

# search for the optimal correspondence between subunits, then calculate CA rmsd
def calculate_rmsd(orientation1, orientation2, homo_chains_list):
    orientation2 = search_subunits_mapping(orientation1, orientation2, homo_chains_list)
    ca_rmsd = calculate_rmsd_orientation(orientation1, orientation2)
    return ca_rmsd


############    MAIN    ############
def main():
    chain_num = int(sys.argv[1])
    csym = int(sys.argv[2])
    nmax = int(sys.argv[3])

    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)

    start_time = time.time()
    # print(f"[INFO] Start C_{csym} symmetric docking", flush=True)

    num_cpus = get_num_cpus(20)

    if csym < chain_num:
        # local symmetry sampling
        comple_csym = chain_num // csym
        mono_file = "A.pdb"
        cn_out_file1 = f"A_c{comple_csym}.out"
        cn_pdb_file1 = f"A{comple_csym}.pdb"
        log_file1 = f"A{comple_csym}.log"
        if os.path.exists(f"{parent_dir}/{cn_out_file1}"):
            shutil.copy(f"{parent_dir}/{cn_out_file1}", current_dir)
            subprocess.run(
                f"{compcn} {cn_out_file1} {cn_pdb_file1} -nmax {nmax} -complex > /dev/null &&"
                f"python {splitmodels_script} {cn_pdb_file1}", shell=True, check=True
            )
        else:
            subprocess.run(
                f"{chdock} {mono_file} {mono_file} -symmetry {comple_csym} -out {cn_out_file1} > {log_file1} &&"
                f"{compcn} {cn_out_file1} {cn_pdb_file1} -nmax {nmax} -complex > /dev/null &&"
                f"python {splitmodels_script} {cn_pdb_file1}", shell=True, check=True
            )

        # global symmetry sampling
        stem = f"A{comple_csym}"
        with ProcessPoolExecutor(max_workers=num_cpus) as executor:
            futures = [executor.submit(run_csym_dock, stem, i, csym, nmax) for i in range(1, nmax+1)]
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    print(f"[ERROR] Failed in global symmetry sampling: {e}")
                    sys.exit(1)
        gc.collect()

        # models clustering
        models = []
        num_models = nmax * nmax
        for i in range(num_models):
            model_file = f"{current_dir}/model_{i + 1}.pdb"
            if Path(model_file).is_file():
                models.append({"idx": i, "file": model_file})
            else:
                raise FileNotFoundError(f"[ERROR] Missing model file: {model_file}.")

        models = calculate_itscores_all(models, num_cpus)
        models.sort(key=lambda x: x["it_score"])

        homo_chain_ids = []
        model1_file = models[0]["file"]
        model1_struct = pdb_parser.get_structure("", model1_file)[0]
        for chain in model1_struct.get_chains():
            homo_chain_ids.append(chain.id)
        for model in models:
            model_struct = pdb_parser.get_structure("", model["file"])[0]
            orientation = {}
            for chain_id in homo_chain_ids:
                chain = model_struct[chain_id]
                orientation[chain_id] = chain
            model["orientation"] = orientation

        clstr_cn_models = []
        rmsd_threshold = 5.0
        for model1 in models:
            added_to_clstr = False
            for model2 in clstr_cn_models:
                ca_rmsd = calculate_rmsd(
                    model1["orientation"], model2["orientation"], homo_chain_ids
                )
                if ca_rmsd < rmsd_threshold:
                    added_to_clstr = True
                    break
            if not added_to_clstr:
                clstr_cn_models.append(model1)
            if len(clstr_cn_models) >= nmax:
                break

        with open(f"{current_dir}/aa2.log", 'w+') as fout:
            for i, cn_model in enumerate(clstr_cn_models):
                file1 = cn_model["file"]
                file2 = f"{current_dir}/cn_{i + 1}.pdb"
                os.rename(file1, file2)
                it_score = cn_model["it_score"]
                print(f"C{csym} model {i + 1} IT-score: {it_score}", file=fout)

        pattern = os.path.join(current_dir, "model_*")
        model_files = glob.glob(pattern)
        for file in model_files:
            try:
                os.remove(file)
            except Exception as e:
                print(f"[ERROR] Failed in deleting {file}: {e}")
                sys.exit(1)

        end_time = time.time()
        print(f"\nC{csym} symmetric docking cost {end_time - start_time}", flush=True)

    elif csym == chain_num:
        # global symmetry sampling
        mono_file = "A.pdb"
        cn_out_file = f"A_c{csym}.out"
        cn_pdb_file = f"A{csym}.pdb"
        log_file = f"A{csym}.log"
        subprocess.run(
            f"{chdock} {mono_file} {mono_file} -symmetry {csym} -out {cn_out_file} > {log_file} &&"
            f"{compcn} {cn_out_file} {cn_pdb_file} -nmax {nmax} -complex > /dev/null &&"
            f"python {splitmodels_script} {cn_pdb_file}", shell=True, check=True
        )

        # models clustering
        stem = f"A{csym}"
        models = []
        for i in range(nmax):
            model_file = f"{current_dir}/{stem}_{i + 1}.pdb"
            if Path(model_file).is_file():
                models.append({"idx": i, "file": model_file})
            else:
                raise FileNotFoundError(f"[ERROR] Missing model file: {model_file}.")
        models = calculate_itscores_all(models, num_cpus)
        models.sort(key=lambda x: x["it_score"])

        end_time = time.time()
        with open(f"{current_dir}/aa2.log", 'w+') as fout:
            for i, cn_model in enumerate(models):
                file1 = cn_model["file"]
                file2 = f"{current_dir}/cn_{i + 1}.pdb"
                os.rename(file1, file2)
                it_score = cn_model["it_score"]
                print(f"Model {i + 1}, IT-score: {it_score}", file=fout)
            print(f"\nC{csym} symmetric docking cost {end_time - start_time}",
                  file=fout, flush=True)

        print(f"\nC{csym} symmetric docking cost {end_time - start_time}", flush=True)

    else:
        raise ValueError(f"[ERROR] Incorrect cyclic symmetry order.")


if __name__ == "__main__":
    main()
