# define the GA process for optimizing medium- / low-quality structures

import math
import random
import gc
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import xl_config as config
import xl_assemble_utils as utils
from xl_assembler import XLink_Assembler


class XLink_GeneticAlgorithm:
    rmsd_threshold = None
    num_cpus = None
    generations = None
    model_num = None
    output = None
    md_steps = None
    md_engine = None
    if_score_threshold1 = None
    if_score_threshold2 = None

    def __init__(self, population, confidence, q_level, fout,
                 population_size=100, mutation_rate=0.1, xl_threshold=0.8):
        self.population = population
        self.confidence = confidence
        self.level = q_level
        self.fout = fout
        current_size = len(self.population)
        if current_size < population_size:
            print(f"\nWarning: adjust population size to {current_size} due to insufficient structures",
                  file=self.fout, flush=True)
        self.population_size = min(current_size, population_size)
        self.mutation_rate = mutation_rate
        self.xl_threshold = xl_threshold
        self.all_individuals = []

    # def the evolve process in GA
    def evolve(self, generation):
        new_individuals = []
        total_log_lines = []
        total_tasks = 2 * self.population_size
        random.seed(config.RANDOM_SEED + generation)
        with ProcessPoolExecutor(max_workers=self.num_cpus) as executor:
            futures = []
            for i in range(total_tasks):
                parent1, parent2 = random.sample(self.population, 2)
                seed = generation * self.population_size + i
                futures.append(
                    executor.submit(
                        individual_generation_task_xl, parent1, parent2,
                        self.mutation_rate, seed, self.xl_threshold
                    )
                )
            for future in as_completed(futures):
                new_results, new_log_lines = future.result()
                new_individuals.extend(new_results)
                total_log_lines.extend(new_log_lines)
        gc.collect()
        for i in total_log_lines:
            print(i, file=self.fout, flush=True)
        return new_individuals

    # check whether the iteration should be terminated after each round
    def check_iteration_stop(self, new_individuals):
        all_results = self.all_individuals + new_individuals
        confidence, xl_scale_factor = utils.check_population_confidence(
            all_results, 20, self.rmsd_threshold
        )
        confidence_threshold1 = self.if_score_threshold1 * xl_scale_factor
        confidence_threshold2 = self.if_score_threshold2 * xl_scale_factor
        if confidence >= confidence_threshold1:
            return confidence, xl_scale_factor, True
        elif confidence >= confidence_threshold2 and self.level == 2:
            return confidence, xl_scale_factor, True
        else:
            return confidence, xl_scale_factor, False

    # add the newly generated structures to the existing population and cluster
    def cluster_results(self, new_individuals):
        score_threshold = self.all_individuals[self.population_size - 1]["if_score"]
        new_individual_num = 0
        population_num = len(self.all_individuals)
        for result1 in new_individuals:
            if result1["if_score"] < score_threshold:
                continue
            if "full_orientation" not in result1:
                full_orientation = utils.create_full_orientation(result1["orientation"])
                result1["full_orientation"] = full_orientation
            added_to_clstr = False
            for i in range(population_num):
                result2 = self.all_individuals[i]
                ca_rmsd = utils.calculate_rmsd(
                    result1["full_orientation"], result2["full_orientation"]
                )
                if ca_rmsd < self.rmsd_threshold:
                    added_to_clstr = True
                    if_score1, if_score2 = result1["if_score"], result2["if_score"]
                    xl_rate1, xl_rate2 = result1["xl_satisfaction"], result2["xl_satisfaction"]
                    if if_score1 > if_score2:
                        self.all_individuals[i] = result1
                    elif if_score1 == if_score2 and xl_rate1 > xl_rate2:
                        self.all_individuals[i] = result1
                    break
            if not added_to_clstr:
                self.all_individuals.append(result1)
                new_individual_num += 1
        print(f"New individual generated in this generation: {new_individual_num}",
              file=self.fout, flush=True)
        self.all_individuals.sort(key=lambda x: x["xl_satisfaction"], reverse=True)
        self.all_individuals.sort(key=lambda x: x["if_score"], reverse=True)
        self.all_individuals = self.all_individuals[:min(2 * self.population_size, len(self.all_individuals))]
        new_population = self.all_individuals[:self.population_size]
        return new_population, new_individual_num

    # define the complete GA iteration process
    def run(self):
        self.all_individuals.extend(self.population)
        no_new_individual = 0
        no_confidence_improvement = 0
        generation = 0
        for k in range(self.generations):
            generation += 1
            stime = datetime.now()
            print("\n**********************************************************", file=self.fout)
            print(f"GA sampling, generation {generation} / {self.generations}\n", file=self.fout, flush=True)
            print(f"\nGA sampling, generation {generation} / {self.generations}", flush=True)
            new_individuals = self.evolve(generation)
            confidence, xl_scale_factor, iteration_stop = self.check_iteration_stop(new_individuals)
            self.population, new_individual_num = self.cluster_results(new_individuals)
            etime = datetime.now()
            print(f"\nConfidence score in generation {generation}: {confidence}", file=self.fout)
            print(f"XL_scale_factor in generation {generation}: {xl_scale_factor}", file=self.fout)
            print(f"\nGeneration {generation} cost {etime - stime}", file=self.fout, flush=True)
            print(f"Generation {generation} cost {etime - stime}", flush=True)
            print("**********************************************************", file=self.fout)
            if iteration_stop:
                break
            if new_individual_num == 0:
                no_new_individual += 1
            else:
                no_new_individual = 0
            if abs(confidence - self.confidence) <= 1e-5:
                no_confidence_improvement += 1
            else:
                no_confidence_improvement = 0
            self.confidence = confidence
            if no_new_individual >= 3 or no_confidence_improvement >= 3:
                print("No new individuals for 3 consecutive generations, stop iteration",
                      file=self.fout, flush=True)
                break
        utils.add_results_itscore_xl(
            self.population, self.num_cpus, self.md_steps, self.md_engine, self.fout
        )
        utils.print_results_itscore(
            self.population, self.output, self.model_num, generation, self.fout
        )


# define the one-side Crossover process:
# randomly remove up to half of the subunits from each parent,
# and reconstructs the complex by reapplying transformations from parents
def cross_attempt_xl(parent, parent_paths, assembler: XLink_Assembler):
    initial_removed_chain_num = math.floor(config.chain_num / 2)
    max_rounds = min(initial_removed_chain_num, 5)
    child, log_line = None, None
    for nround in range(max_rounds):
        removed_chain_num = initial_removed_chain_num - nround
        if removed_chain_num <= 0:
            break
        max_attempts = min(math.comb(config.chain_num, removed_chain_num), 10)
        removed_chain_ids = random.sample(
            list(combinations(config.chain_ids, removed_chain_num)), max_attempts
        )
        removed_attempt = 0
        while removed_attempt < max_attempts:
            removed_chain_id = removed_chain_ids[removed_attempt]
            removed_attempt += 1
            old_path = [
                interface for interface in parent["path"]
                if not (interface["chain1"] in removed_chain_id or
                        interface["chain2"] in removed_chain_id)
            ]
            checked_interfaces = assembler.check_old_path(old_path)
            if not checked_interfaces:
                continue
            assemblied_all, xl_excluded_ifs = assembler.assembly_in_order(parent_paths, checked_interfaces)
            if assemblied_all:
                assembler.orientation = {
                    chain_id: assembler.orientation[chain_id] for chain_id in config.chain_ids
                }
                chain_id = config.chain_ids[0]
                if_score = assembler.calculate_interface_score(chain_id)
                xl_satisfaction = utils.calculate_xl_satisfaction(config.xl_total_data, assembler.xl_ok_data)
                new_if_score = if_score * xl_satisfaction
                child = {
                    "if_score": new_if_score, "path": assembler.path,
                    "orientation": assembler.orientation, "xl_satisfaction": xl_satisfaction,
                    "original_if_score": if_score
                }
                break
        if child:
            break
    if not child:
        child = parent
        log_line = "Too many attempts, use parent structure"
    return child, log_line


# define the Crossover operation for GA
def ga_crossover_xl(parent1, parent2, assembler: XLink_Assembler):
    parent_paths = [
        dict(i) for i in set(tuple(p.items()) for p in parent1["path"] + parent2["path"])
    ]
    parent_paths = utils.create_path_df(config.interfaces, parent_paths)
    parent_paths = utils.sort_dataframe(
        parent_paths, "scaled_score2", "score", "chain1", "chain2"
    )
    child1, log_line1 = cross_attempt_xl(parent1, parent_paths, assembler)
    child2, log_line2 = cross_attempt_xl(parent2, parent_paths, assembler)
    return child1, child2, log_line1, log_line2


# define the Mutate operation for GA
def ga_mutate_xl(individual, mutation_rate, assembler: XLink_Assembler):
    if random.random() > mutation_rate:
        return individual, None, None
    removed_chain_id = random.choice(config.chain_ids)
    old_path = [
        interface for interface in individual["path"]
        if removed_chain_id not in (interface["chain1"], interface["chain2"])
    ]
    old_path = assembler.check_old_path(old_path)
    rest_interfaces = [
        interface for _, interface in config.interfaces.iterrows()
        if not
        (interface["chain1"] in assembler.groups and
         interface["chain2"] in assembler.groups[interface["chain1"]])
    ]
    fill_attempt = 0
    log_line = None
    while fill_attempt < 200:
        fill_attempt += 1
        interface = random.choice(rest_interfaces)
        chain1, chain2 = interface["chain1"], interface["chain2"]
        if chain1 not in assembler.groups and chain2 not in assembler.groups:
            edge_available, xl_ok = assembler.clash_new(chain1, chain2, interface)
            if not edge_available or not xl_ok:
                continue
        elif chain1 not in assembler.groups:
            edge_available, xl_ok = assembler.clash_append(chain1, chain2, interface, 0.02)
            if not edge_available or not xl_ok:
                continue
        elif chain2 not in assembler.groups:
            edge_available, xl_ok = assembler.clash_append(chain2, chain1, interface, 0.02)
            if not edge_available or not xl_ok:
                continue
        elif chain1 not in assembler.groups[chain2]:
            edge_available, xl_ok = assembler.clash_merge(chain1, chain2, interface)
            if not edge_available or not xl_ok:
                continue
        else:
            continue
        if len(assembler.groups[chain1]) == config.chain_num:
            assembler.assemblied_all = True
            break
    if assembler.assemblied_all:
        assembler.orientation = {
            chain_id: assembler.orientation[chain_id] for chain_id in config.chain_ids
        }
        chain_id = config.chain_ids[0]
        if_score = assembler.calculate_interface_score(chain_id)
        xl_satisfaction = utils.calculate_xl_satisfaction(config.xl_total_data, assembler.xl_ok_data)
        new_if_score = if_score * xl_satisfaction
        new_individual = {
            "if_score": new_if_score, "path": assembler.path, "orientation": assembler.orientation,
            "xl_satisfaction": xl_satisfaction, "original_if_score": if_score
        }
        return individual, new_individual, log_line
    else:
        log_line = "Too many attempts, use original individual"
        return individual, None, log_line


# define the generation process of the new structures, by Crossover and Mutate
def individual_generation_task_xl(parent1, parent2, mutation_rate, seed, xl_threshold):
    random.seed(seed)
    single_assembler = XLink_Assembler(xl_threshold)
    log_lines = []
    child1, child2, log_line1, log_line2 = ga_crossover_xl(parent1, parent2, single_assembler)
    log_lines.extend([log_line1, log_line2])
    child1, new_individual1, log_line1 = ga_mutate_xl(child1, mutation_rate, single_assembler)
    child2, new_individual2, log_line2 = ga_mutate_xl(child2, mutation_rate, single_assembler)
    log_lines.extend([log_line1, log_line2])
    del single_assembler
    gc.collect()
    new_results = [child1, child2, new_individual1, new_individual2]
    filtered_results = [i for i in new_results if i]
    filtered_log_lines = [i for i in log_lines if i]
    return filtered_results, filtered_log_lines
