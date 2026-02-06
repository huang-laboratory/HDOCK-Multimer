import gc
import sys
import os
import glob
import time
import subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, "..", ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))

sys.path.append(BASE_PATH)
scorecom = os.path.join(TOOL_PATH, "scorecom.sh")
dhdock = os.path.join(TOOL_PATH, "dhdock.sh")
compdn = os.path.join(TOOL_PATH, "compdn")
splitmodels_script = os.path.join(TOOL_PATH, "splitmodels.py")


############    Utility Functions    ############
# select appropriate number of parallel CPU cores
def get_num_cpus(max_workers: int):
    cpu_cap = int(os.getenv("NUM_CPUS", max_workers))
    cpu_count = os.cpu_count()
    return max(1, min(cpu_count, cpu_cap))

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


############    MAIN    ############
def main():
    dsym = int(sys.argv[1])
    nmax = int(sys.argv[2])

    current_dir = os.getcwd()

    start_time = time.time()
    # print(f"[INFO] Start D_{dsym} symmetric docking", flush=True)

    num_cpus = get_num_cpus(20)

    # Dn symmetry sampling
    mono_file = "A.pdb"
    dn_out_file = f"A_d{dsym}.out"
    dn_pdb_file = f"A_d{dsym}.pdb"
    log_file = f"A_d{dsym}.log"
    subprocess.run(
        f"{dhdock} {mono_file} -dn {dsym} -out {dn_out_file} -ncmer {nmax} > {log_file} &&"
        f"{compdn} {dn_out_file} {dn_pdb_file} -nmax {nmax} -complex > /dev/null &&"
        f"python {splitmodels_script} {dn_pdb_file}", shell=True, check=True
    )

    # models clustering
    stem = f"A_d{dsym}"
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
        for i, dn_model in enumerate(models):
            file1 = dn_model["file"]
            file2 = f"{current_dir}/dn_{i + 1}.pdb"
            os.rename(file1, file2)
            it_score = dn_model["it_score"]
            print(f"Model {i + 1}, IT-score: {it_score}", file=fout)
        print(f"\nD{dsym} symmetric docking cost {end_time - start_time}",
              file=fout, flush=True)

    print(f"\nD{dsym} symmetric docking cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    main()
