import os
import shutil
from argparse import ArgumentParser
import time
import config
import assemble_utils as utils
import assembler


def main(fout):
    # parse input arguments 
    arg_parser = ArgumentParser()
    arg_parser.add_argument("-d", "--pdbdir", type=str, required=True)
    arg_parser.add_argument("-s", "--stoi", type=str, required=True)
    arg_parser.add_argument("-n", "--nmax", type=int, default=10)
    arg_parser.add_argument("-nr", "--nround", type=int, default=500)
    arg_parser.add_argument("-r", "--rmsd", type=float, default=5.0)
    arg_parser.add_argument("-p", "--population", type=int, default=100)
    arg_parser.add_argument("-w", "--workers", type=int, default=20)
    arg_parser.add_argument("-g", "--generations", type=int, default=50)
    arg_parser.add_argument("-m", "--md", type=int, default=500)
    arg_parser.add_argument("-e", "--md_engine", type=str,
                            choices=["openmm", "amber"], default="openmm")

    args = arg_parser.parse_args()
    subcomplex_dir = args.pdbdir
    file_stoi = args.stoi
    model_num = args.nmax
    sampling_round = args.nround
    rmsd = args.rmsd
    population_size = args.population
    max_workers = args.workers
    num_cpus = config.get_num_cpus(max_workers)
    generations = args.generations
    md_steps = args.md
    md_engine = args.md_engine

    start_time = time.time()
    # load data
    config.load_struct_data(subcomplex_dir, file_stoi, num_cpus, fout)

    # start initail sampling
    assemblied_success, connectivity, initial_results = assembler.generate_initial_population(
        num_cpus, sampling_round, fout
    )
    if assemblied_success:
        confidence = utils.check_population_confidence(initial_results, 20, rmsd)
        print(f"\nConfidence score in generation 0: {confidence}", file=fout, flush=True)
        initial_results = utils.cluster_results(initial_results, population_size, rmsd, 1)
        num_models = min(len(initial_results), model_num)
        for i in range(num_models):
            confidence_score = initial_results[i]["if_score"]
            print(f"Confidence: initial model {i + 1} confidence_score: {confidence_score}",
                  file=fout, flush=True)
        initial_results = utils.add_results_itscore(
            initial_results, num_cpus, md_steps, md_engine, fout
        )
        initial_results.sort(key=lambda x: x["it_score"])
        for i in range(num_models):
            it_score = initial_results[i]["it_score"]
            normalized_it_score = it_score / config.total_length
            print(f"Normalized_ITscore: initial model {i + 1} normalized it_score: {normalized_it_score}",
                  file=fout, flush=True)

    elif connectivity:
        print(f"\nAssembly failed, print optimal subcomplexes for subsequent processing",
              file=fout, flush=True)
        optimal_result = utils.select_optimal_subcomplex(initial_results)
        utils.print_groups(optimal_result, connectivity, fout)

    else:
        print(f"\nAssembly failed, subunits connectivity is not satisfied",
              file=fout, flush=True)
        print(f"Print optimal subcomplexes for subsequent processing",
              file=fout, flush=True)
        optimal_result = utils.select_optimal_subcomplex(initial_results)
        utils.print_groups(optimal_result, connectivity, fout)

    if os.path.isdir(config.temp_dir):
        shutil.rmtree(config.temp_dir)

    end_time = time.time()
    print(f"\nAssembly process cost {end_time - start_time}", file=fout, flush=True)
    print(f"\nAssembly sampling process cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    with open(f"{config.current_dir}/aa2.log", 'w+') as fout:
        main(fout)
