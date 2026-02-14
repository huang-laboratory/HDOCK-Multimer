import os
import shutil
from argparse import ArgumentParser
import time
import xl_config as config
import xl_assemble_utils as utils
import xl_assembler as assembler
from xl_geneticalgorithm import XLink_GeneticAlgorithm


####################        MAIN         ####################
def main(fout):
    # parse input parameters
    arg_parser = ArgumentParser()
    arg_parser.add_argument("-o", "--output", type=str, required=True)
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
    arg_parser.add_argument("-s1", "--score1", type=float, default=85.0)
    arg_parser.add_argument("-s2", "--score2", type=float, default=80.0)
    # added for parsing crosslinking arguments
    arg_parser.add_argument("-x", "--xlinks", type=str)
    arg_parser.add_argument("-xt1", "--xl_threshold1", type=float, default=0.7)
    arg_parser.add_argument("-xt2", "--xl_threshold2", type=float, default=0.8)
    arg_parser.add_argument("-xo", "--xl_offset", type=int, default=20)
    arg_parser.add_argument("-xf", "--xl_self", type=bool, default=False)

    args = arg_parser.parse_args()
    subcomplex_dir = args.pdbdir
    file_output = args.output
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
    if_score_threshold1 = args.score1
    if_score_threshold2 = args.score2

    file_crosslinks = args.xlinks
    xl_threshold1 = args.xl_threshold1
    xl_threshold2 = args.xl_threshold2
    xl_max_offset = args.xl_offset
    allow_self_linked = args.xl_self

    start_time = time.time()

    # load data
    config.load_struct_data(subcomplex_dir, file_stoi, num_cpus, fout)
    config.load_crosslink_data(file_crosslinks, xl_max_offset, allow_self_linked)

    XLink_GeneticAlgorithm.rmsd_threshold = rmsd
    XLink_GeneticAlgorithm.num_cpus = num_cpus
    XLink_GeneticAlgorithm.generations = generations
    XLink_GeneticAlgorithm.model_num = model_num
    XLink_GeneticAlgorithm.output = file_output
    XLink_GeneticAlgorithm.md_steps = md_steps
    XLink_GeneticAlgorithm.md_engine = md_engine
    XLink_GeneticAlgorithm.if_score_threshold1 = if_score_threshold1
    XLink_GeneticAlgorithm.if_score_threshold2 = if_score_threshold2

    # start assembly process
    assemblied_success, connectivity, initial_results = assembler.generate_initial_population_xl(
        num_cpus, sampling_round, fout, xl_threshold1
    )
    if assemblied_success:
        confidence, xl_scale_factor = utils.check_population_confidence(initial_results, 20, rmsd)
        print(f"\nConfidence score in generation 0: {confidence}", file=fout, flush=True)
        print(f"XL_scale_factor in generation 0: {xl_scale_factor}\n", file=fout, flush=True)
        confidence_threshold1 = if_score_threshold1 * xl_scale_factor
        confidence_threshold2 = if_score_threshold2 * xl_scale_factor
        early_stop = False
        initial_results = utils.cluster_results(initial_results, population_size, rmsd, 1)
        if len(initial_results) <= model_num:
            early_stop = True
        if confidence >= confidence_threshold1:
            q_level = 0
        elif confidence >= confidence_threshold2:
            q_level = 1
        else:
            q_level = 2
        if q_level == 0 or early_stop:
            initial_results = utils.add_results_itscore_xl(
                initial_results, num_cpus, md_steps, md_engine, fout
            )
            utils.print_results_itscore(
                initial_results, file_output, model_num, 0, fout
            )
        else:
            print(f"\nInitial population quality level: {q_level}", file=fout, flush=True)
            print(f"\nTry to optimize population via GA iteration", file=fout, flush=True)
            print(f"\nInitial population quality level: {q_level}", flush=True)
            print(f"\nTry to optimize population via GA iteration", flush=True)
            ga_generator = XLink_GeneticAlgorithm(
                initial_results, confidence, q_level, fout, population_size, 0.1, xl_threshold2
            )
            ga_generator.run()

    elif connectivity:
        print(f"\nAssembly failed, print optimal subcomplexes for multi-body docking",
              file=fout, flush=True)
        print(f"\nAssembly failed, print optimal subcomplexes for multi-body docking",
              flush=True)
        optimal_result = utils.select_optimal_subcomplex(initial_results)
        utils.print_groups(optimal_result, connectivity, fout)

    else:
        print("\n[WARN] Assembly failed, subunits connectivity is not satisfied",
              file=fout, flush=True)
        print("[WARN] Print optimal subcomplexes for subsequent processing",
              file=fout, flush=True)
        print("\n[WARN] Assembly failed, subunits connectivity is not satisfied", flush=True)
        optimal_result = utils.select_optimal_subcomplex(initial_results)
        utils.print_groups(optimal_result, connectivity, fout)

    if os.path.isdir(config.temp_dir):
        shutil.rmtree(config.temp_dir)

    end_time = time.time()
    print(f"\nAssembly process cost {end_time - start_time}", file=fout, flush=True)
    print(f"\nAssembly process cost {end_time - start_time}", flush=True)


if __name__ == "__main__":
    with open(f"{config.current_dir}/aa2.log", 'w+') as fout:
        main(fout)

