# define Assembler and the random sampling process

import numpy as np
import gc
import copy
from scipy.spatial.distance import cdist
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import xl_config as config
import xl_assemble_utils as utils


steric_clash_dist = config.steric_clash_dist
BACKBONE = config.BACKBONE


class XLink_Assembler:
    def __init__(self, xl_threshold):
        self.orientation = {}
        self.groups = {}
        self.transforms = {}
        self.assemblied_all = False
        self.path = []
        self.sum_wt = {}
        self.sum_wtst1 = {}
        self.sum_wtst2 = {}
        # added for crosslinking restraints
        self.xl_threshold = xl_threshold
        self.xl_ok_data = {}

    def reset(self):
        self.orientation = {}
        self.groups = {}
        self.transforms = {}
        self.assemblied_all = False
        self.path = []
        self.sum_wt = {chain_id: 0.0 for chain_id in config.chain_ids}
        self.sum_wtst1 = {chain_id: 0.0 for chain_id in config.chain_ids}
        self.sum_wtst2 = {chain_id: 0.0 for chain_id in config.chain_ids}
        # added for crosslinking restraints
        self.xl_ok_data = {}

    def count_group_xls(self, chains, total_xls: dict, satisfied_xls: dict):
        pairs = [
            config.generate_chain_pair(c1, c2) for c1, c2 in combinations(chains, 2)
        ]
        for pair in pairs:
            if pair in config.pair_restraint_ids:
                for xl_id in config.pair_restraint_ids[pair]:
                    restraint = config.xl_restraints[xl_id]["xl_pairs"][pair]
                    w1, w2 = restraint["w1"], restraint["w2"]
                    xl_data = {"pair": pair, "w1": w1, "w2": w2}
                    xl_existed = xl_id in self.xl_ok_data
                    utils.upsert_xl_max_w2(total_xls, xl_id, xl_data)
                    if xl_existed:
                        utils.upsert_xl_max_w2(satisfied_xls, xl_id, xl_data)
        return total_xls, satisfied_xls

    # create: assemble two unpaired subunits
    def clash_new(self, chain1, chain2, interface, ratio: float = 0.02):
        source = interface["source"]
        mono_chain1 = config.mono_structs[chain1].copy()
        mono_chain2 = config.mono_structs[chain2].copy()
        rot1, trans1 = config.transformations[chain1][source]
        rot2, trans2 = config.transformations[chain2][source]
        for atom in mono_chain1.get_atoms():
            atom.transform(rot1, trans1)
        for atom in mono_chain2.get_atoms():
            atom.transform(rot2, trans2)

        chain1_atoms = np.array([
            atom.get_coord() for atom in mono_chain1.get_atoms()
            if atom.get_name() in BACKBONE and atom.get_bfactor() > 80
        ])
        chain2_atoms = np.array([
            atom.get_coord() for atom in mono_chain2.get_atoms()
            if atom.get_name() in BACKBONE and atom.get_bfactor() > 80
        ])
        steric_clash = False
        subunits_connectivity = True
        if len(chain1_atoms) > 0 and len(chain2_atoms) > 0:
            chain1_atoms_num = sum(1 for atom in mono_chain1.get_atoms() if atom.get_name() in BACKBONE)
            chain2_atoms_num = sum(1 for atom in mono_chain2.get_atoms() if atom.get_name() in BACKBONE)
            clashes = np.argwhere(cdist(chain1_atoms, chain2_atoms) < steric_clash_dist)
            clash1, clash2 = set(clashes[:, 0]), set(clashes[:, 1])
            if len(clash1) > chain1_atoms_num * ratio or len(clash2) > chain2_atoms_num * ratio:
                steric_clash = True
        if chain1 in config.subunits_id_mapping:
            full_chain_id = config.subunits_id_mapping[chain1]
            order_map = config.subunits_order[full_chain_id]
            adjacent = chain2 in order_map and abs(order_map[chain1] - order_map[chain2]) == 1
            if adjacent:
                subunits_connectivity = utils.check_subunits_connectivity(
                    chain1, chain2, mono_chain1, mono_chain2
                )
        edge_available = subunits_connectivity and not steric_clash
        if not edge_available:
            return edge_available, False
        else:
            total_xls, satisfied_xls = {}, {}
            newly_ok_xls = {}
            pair = config.generate_chain_pair(chain1, chain2)
            if pair in config.pair_restraint_ids:
                for xl_id in config.pair_restraint_ids[pair]:
                    restraint = config.xl_restraints[xl_id]["xl_pairs"][pair]
                    w1, w2 = restraint["w1"], restraint["w2"]
                    xl_data = {"pair": pair, "w1": w1, "w2": w2}
                    total_xls[xl_id] = xl_data
                    xl_existed = xl_id in self.xl_ok_data
                    xl_pass = utils.check_xl_for_pair(
                        chain1, chain2, mono_chain1, mono_chain2, restraint
                    )
                    if xl_pass:
                        satisfied_xls[xl_id] = xl_data
                        newly_ok_xls[xl_id] = xl_data
                    elif xl_existed:
                        satisfied_xls[xl_id] = xl_data
            xl_satisfaction = utils.calculate_xl_satisfaction(total_xls, satisfied_xls)
            xl_ok = xl_satisfaction >= self.xl_threshold
            if xl_ok:
                aa_num1 = int(interface["length1"])
                aa_num2 = int(interface["length2"])
                weight = min(aa_num1, aa_num2)
                scaled_score1 = 0.4 * interface["scaled_score1"] + 0.6 * interface["scaled_score2"]
                scaled_score2 = interface["scaled_score2"]
                self.sum_wt[chain1] += weight
                self.sum_wtst1[chain1] += weight * scaled_score1
                self.sum_wtst2[chain1] += weight * scaled_score2
                self.sum_wt[chain2] += weight
                self.sum_wtst1[chain2] += weight * scaled_score1
                self.sum_wtst2[chain2] += weight * scaled_score2
                self.orientation[chain1] = mono_chain1
                self.transforms[chain1] = [(rot1, trans1)]
                self.orientation[chain2] = mono_chain2
                self.transforms[chain2] = [(rot2, trans2)]
                new_complex = [chain1, chain2]
                self.groups[chain1] = new_complex
                self.groups[chain2] = new_complex
                self.path.append(interface)

                for xl_id, xl_data in newly_ok_xls.items():
                    utils.upsert_xl_max_w2(self.xl_ok_data, xl_id, xl_data)
            return edge_available, xl_ok

    # append: incorporate an unpaired subunit in an existed subcomplex
    def clash_append(self, chain1, chain2, interface, ratio: float = 0.05):
        source = interface["source"]
        mono_chain1 = config.mono_structs[chain1].copy()
        rot1, trans1 = config.transformations[chain1][source]
        rot2, trans2 = config.transformations[chain2][source]
        for atom in mono_chain1.get_atoms():
            atom.transform(rot1, trans1)
        inverse_transforms = utils.reverse_transform_list(self.transforms[chain2])
        inverse_transforms.append((rot2, trans2))
        total_rot, total_trans = utils.combine_transform_list(inverse_transforms)
        for chain in self.groups[chain2]:
            for atom in self.orientation[chain].get_atoms():
                atom.transform(total_rot, total_trans)

        chain1_atoms_num = len(list(mono_chain1.get_atoms()))
        chain1_atoms = np.array([
            atom.get_coord() for atom in mono_chain1.get_atoms()
            if atom.get_name() in BACKBONE and atom.get_bfactor() > 80
        ])
        chain2_atoms = np.array([
            atom.get_coord()
            for chain in self.groups[chain2]
            for atom in self.orientation[chain].get_atoms()
            if atom.get_name() in BACKBONE and atom.get_bfactor() > 80
        ])
        steric_clash = False
        subunits_connectivity = True
        if len(chain1_atoms) > 0 and len(chain2_atoms) > 0:
            clashes = np.argwhere(cdist(chain1_atoms, chain2_atoms) < steric_clash_dist)
            clash_chain1 = set(clashes[:, 0])
            if len(clash_chain1) > chain1_atoms_num * ratio:
                steric_clash = True
        if chain1 in config.subunits_id_mapping:
            full_chain_id = config.subunits_id_mapping[chain1]
            order_map = config.subunits_order[full_chain_id]
            for chain in self.groups[chain2]:
                adjacent = chain in order_map and abs(order_map[chain1] - order_map[chain]) == 1
                if adjacent:
                    subunits_connectivity = utils.check_subunits_connectivity(
                        chain1, chain, mono_chain1, self.orientation[chain]
                    )
                    if not subunits_connectivity:
                        break

        edge_available = subunits_connectivity and not steric_clash
        if not edge_available:
            inverse_rot, inverse_trans = utils.reverse_transform(total_rot, total_trans)
            for chain in self.groups[chain2]:
                for atom in self.orientation[chain].get_atoms():
                    atom.transform(inverse_rot, inverse_trans)
            return edge_available, False
        else:
            total_xls, satisfied_xls = {}, {}
            total_xls, satisfied_xls = self.count_group_xls(self.groups[chain2], total_xls, satisfied_xls)
            newly_ok_xls = {}
            for chain in self.groups[chain2]:
                pair = config.generate_chain_pair(chain1, chain)
                if pair in config.pair_restraint_ids:
                    for xl_id in config.pair_restraint_ids[pair]:
                        restraint = config.xl_restraints[xl_id]["xl_pairs"][pair]
                        w1, w2 = restraint["w1"], restraint["w2"]
                        xl_data = {"pair": pair, "w1": w1, "w2": w2}
                        utils.upsert_xl_max_w2(total_xls, xl_id, xl_data)
                        xl_existed = xl_id in self.xl_ok_data
                        xl_pass = utils.check_xl_for_pair(
                            chain1, chain, mono_chain1, self.orientation[chain], restraint
                        )
                        if xl_pass:
                            utils.upsert_xl_max_w2(satisfied_xls, xl_id, xl_data)
                            utils.upsert_xl_max_w2(newly_ok_xls, xl_id, xl_data)
                        elif xl_existed:
                            utils.upsert_xl_max_w2(satisfied_xls, xl_id, xl_data)
            xl_satisfaction = utils.calculate_xl_satisfaction(total_xls, satisfied_xls)
            xl_ok = xl_satisfaction >= self.xl_threshold
            if not xl_ok:
                inverse_rot, inverse_trans = utils.reverse_transform(total_rot, total_trans)
                for chain in self.groups[chain2]:
                    for atom in self.orientation[chain].get_atoms():
                        atom.transform(inverse_rot, inverse_trans)
            else:
                aa_num1 = sum(1 for residue in mono_chain1.get_residues())
                aa_num2 = sum(
                    sum(1 for residue in self.orientation[chain].get_residues())
                    for chain in self.groups[chain2]
                )
                weight = min(aa_num1, aa_num2)
                scaled_score1 = 0.4 * interface["scaled_score1"] + 0.6 * interface["scaled_score2"]
                scaled_score2 = interface["scaled_score2"]
                self.orientation[chain1] = mono_chain1
                self.transforms[chain1] = [(rot1, trans1)]
                for chain in self.groups[chain2]:
                    self.transforms[chain].append((total_rot, total_trans))
                self.groups[chain2].append(chain1)
                self.groups[chain1] = self.groups[chain2]
                sum_wt = self.sum_wt[chain2] + weight
                sum_wtst1 = self.sum_wtst1[chain2] + weight * scaled_score1
                sum_wtst2 = self.sum_wtst2[chain2] + weight * scaled_score2
                for i in self.groups[chain2]:
                    self.sum_wt[i] = sum_wt
                    self.sum_wtst1[i] = sum_wtst1
                    self.sum_wtst2[i] = sum_wtst2
                self.path.append(interface)

                for xl_id, xl_data in newly_ok_xls.items():
                    utils.upsert_xl_max_w2(self.xl_ok_data, xl_id, xl_data)
            return edge_available, xl_ok

    # merge: combine two different subcomplexes into a larger complex
    def clash_merge(self, chain1, chain2, interface, ratio: float = 0.02):
        source = interface["source"]
        rot1, trans1 = config.transformations[chain1][source]
        rot2, trans2 = config.transformations[chain2][source]
        inverse_transforms1 = utils.reverse_transform_list(self.transforms[chain1])
        inverse_transforms2 = utils.reverse_transform_list(self.transforms[chain2])
        inverse_transforms1.append((rot1, trans1))
        inverse_transforms2.append((rot2, trans2))
        total_rot1, total_trans1 = utils.combine_transform_list(inverse_transforms1)
        total_rot2, total_trans2 = utils.combine_transform_list(inverse_transforms2)
        for chain in self.groups[chain1]:
            for atom in self.orientation[chain].get_atoms():
                atom.transform(total_rot1, total_trans1)
        for chain in self.groups[chain2]:
            for atom in self.orientation[chain].get_atoms():
                atom.transform(total_rot2, total_trans2)

        chain1_atoms = np.array([
            atom.get_coord()
            for chain in self.groups[chain1]
            for atom in self.orientation[chain].get_atoms()
            if atom.get_name() in BACKBONE and atom.get_bfactor() > 80
        ])
        chain2_atoms = np.array([
            atom.get_coord()
            for chain in self.groups[chain2]
            for atom in self.orientation[chain].get_atoms()
            if atom.get_name() in BACKBONE and atom.get_bfactor() > 80
        ])
        steric_clash = False
        subunits_connectivity = True
        if len(chain1_atoms) > 0 and len(chain2_atoms) > 0:
            chain1_atoms_num = sum(
                sum(1 for atom in self.orientation[chain].get_atoms() if atom.get_name() in BACKBONE)
                for chain in self.groups[chain1]
            )
            chain2_atoms_num = sum(
                sum(1 for atom in self.orientation[chain].get_atoms() if atom.get_name() in BACKBONE)
                for chain in self.groups[chain2]
            )
            clashes = np.argwhere(cdist(chain1_atoms, chain2_atoms) < steric_clash_dist)
            clash1, clash2 = set(clashes[:, 0]), set(clashes[:, 1])
            if chain1_atoms_num < chain2_atoms_num:
                steric_clash = len(clash1) > chain1_atoms_num * ratio
            else:
                steric_clash = len(clash2) > chain2_atoms_num * ratio
        for chain_id1 in self.groups[chain1]:
            if chain_id1 in config.subunits_id_mapping:
                full_chain_id1 = config.subunits_id_mapping[chain_id1]
                order_map1 = config.subunits_order[full_chain_id1]
                for chain_id2 in self.groups[chain2]:
                    adjacent = chain_id2 in order_map1 and abs(order_map1[chain_id1] - order_map1[chain_id2]) == 1
                    if adjacent:
                        subunits_connectivity = utils.check_subunits_connectivity(
                            chain_id1, chain_id2, self.orientation[chain_id1], self.orientation[chain_id2]
                        )
                        if not subunits_connectivity:
                            break
                if not subunits_connectivity:
                    break

        edge_available = subunits_connectivity and not steric_clash
        if not edge_available:
            inverse_rot1, inverse_trans1 = utils.reverse_transform(total_rot1, total_trans1)
            inverse_rot2, inverse_trans2 = utils.reverse_transform(total_rot2, total_trans2)
            for chain in self.groups[chain1]:
                for atom in self.orientation[chain].get_atoms():
                    atom.transform(inverse_rot1, inverse_trans1)
            for chain in self.groups[chain2]:
                for atom in self.orientation[chain].get_atoms():
                    atom.transform(inverse_rot2, inverse_trans2)
            return edge_available, False
        else:
            total_xls, satisfied_xls = {}, {}
            total_xls, satisfied_xls = self.count_group_xls(self.groups[chain1], total_xls, satisfied_xls)
            total_xls, satisfied_xls = self.count_group_xls(self.groups[chain2], total_xls, satisfied_xls)
            newly_ok_xls = {}
            for chain_id1 in self.groups[chain1]:
                for chain_id2 in self.groups[chain2]:
                    pair = config.generate_chain_pair(chain_id1, chain_id2)
                    if pair in config.pair_restraint_ids:
                        for xl_id in config.pair_restraint_ids[pair]:
                            restraint = config.xl_restraints[xl_id]["xl_pairs"][pair]
                            w1, w2 = restraint["w1"], restraint["w2"]
                            xl_data = {"pair": pair, "w1": w1, "w2": w2}
                            utils.upsert_xl_max_w2(total_xls, xl_id, xl_data)
                            xl_existed = xl_id in self.xl_ok_data
                            xl_pass = utils.check_xl_for_pair(
                                chain_id1, chain_id2, self.orientation[chain_id1],
                                self.orientation[chain_id2], restraint
                            )
                            if xl_pass:
                                utils.upsert_xl_max_w2(satisfied_xls, xl_id, xl_data)
                                utils.upsert_xl_max_w2(newly_ok_xls, xl_id, xl_data)
                            elif xl_existed:
                                utils.upsert_xl_max_w2(satisfied_xls, xl_id, xl_data)
            xl_satisfaction = utils.calculate_xl_satisfaction(total_xls, satisfied_xls)
            xl_ok = xl_satisfaction >= self.xl_threshold
            if not xl_ok:
                inverse_rot1, inverse_trans1 = utils.reverse_transform(total_rot1, total_trans1)
                inverse_rot2, inverse_trans2 = utils.reverse_transform(total_rot2, total_trans2)
                for chain in self.groups[chain1]:
                    for atom in self.orientation[chain].get_atoms():
                        atom.transform(inverse_rot1, inverse_trans1)
                for chain in self.groups[chain2]:
                    for atom in self.orientation[chain].get_atoms():
                        atom.transform(inverse_rot2, inverse_trans2)
            else:
                aa_num1 = sum(
                    sum(1 for residue in self.orientation[chain].get_residues())
                    for chain in self.groups[chain1]
                )
                aa_num2 = sum(
                    sum(1 for residue in self.orientation[chain].get_residues())
                    for chain in self.groups[chain2]
                )
                weight = min(aa_num1, aa_num2)
                scaled_score1 = 0.4 * interface["scaled_score1"] + 0.6 * interface["scaled_score2"]
                scaled_score2 = interface["scaled_score2"]
                for chain in self.groups[chain1]:
                    self.transforms[chain].append((total_rot1, total_trans1))
                for chain in self.groups[chain2]:
                    self.transforms[chain].append((total_rot2, total_trans2))
                new_complex = self.groups[chain1] + self.groups[chain2]
                sum_wt = self.sum_wt[chain1] + self.sum_wt[chain2] + weight
                sum_wtst1 = self.sum_wtst1[chain1] + self.sum_wtst1[chain2] + weight * scaled_score1
                sum_wtst2 = self.sum_wtst2[chain1] + self.sum_wtst2[chain2] + weight * scaled_score2
                for i in new_complex:
                    self.groups[i] = new_complex
                    self.sum_wt[i] = sum_wt
                    self.sum_wtst1[i] = sum_wtst1
                    self.sum_wtst2[i] = sum_wtst2
                self.path.append(interface)

                for xl_id, xl_data in newly_ok_xls.items():
                    utils.upsert_xl_max_w2(self.xl_ok_data, xl_id, xl_data)
            return edge_available, xl_ok

    # calculate weighted interface score (confidence score)
    def calculate_interface_score(self, chain):
        if self.sum_wt[chain] != 0.0:
            score = self.sum_wtst2[chain] / self.sum_wt[chain]
        else:
            score = 0.0
        return score

    # not used in this study
    def calculate_snew_score(self, chain):
        if self.sum_wt[chain] != 0.0:
            score = self.sum_wtst1[chain] / self.sum_wt[chain]
        else:
            score = 0.0
        return score

    # traverse the randomly shuffled transformations to assemble full complex
    def assembly_in_order(self, interfaces, removed_ifs=None):
        xl_excluded_ifs = []
        for _, interface in interfaces.iterrows():
            chain1, chain2, source = interface["chain1"], interface["chain2"], interface["source"]
            if removed_ifs is not None:
                if (chain1, chain2, source) in removed_ifs:
                    continue
            if chain1 not in self.groups and chain2 not in self.groups:
                edge_available, xl_ok = self.clash_new(chain1, chain2, interface)
                if not edge_available:
                    continue
                elif not xl_ok:
                    xl_excluded_ifs.append(interface)
                    continue
            elif chain1 not in self.groups:
                edge_available, xl_ok = self.clash_append(chain1, chain2, interface)
                if not edge_available:
                    continue
                elif not xl_ok:
                    xl_excluded_ifs.append(interface)
                    continue
            elif chain2 not in self.groups:
                edge_available, xl_ok = self.clash_append(chain2, chain1, interface)
                if not edge_available:
                    continue
                elif not xl_ok:
                    xl_excluded_ifs.append(interface)
                    continue
            elif chain1 not in self.groups[chain2]:
                edge_available, xl_ok = self.clash_merge(chain1, chain2, interface)
                if not edge_available:
                    continue
                elif not xl_ok:
                    xl_excluded_ifs.append(interface)
                    continue
            else:
                continue
            if len(self.groups[chain1]) == config.chain_num:
                self.assemblied_all = True
                break
        return self.assemblied_all, xl_excluded_ifs

    # assemble along the old path until a conflict occurs
    def check_old_path(self, old_path):
        self.reset()
        checked_ifs = []
        for interface in old_path:
            chain1, chain2, source = interface["chain1"], interface["chain2"], interface["source"]
            if chain1 not in self.groups and chain2 not in self.groups:
                edge_available, xl_ok = self.clash_new(chain1, chain2, interface)
                if not edge_available:
                    checked_ifs.append((chain1, chain2, source))
                    break
                elif not xl_ok:
                    break
            elif chain1 not in self.groups:
                edge_available, xl_ok = self.clash_append(chain1, chain2, interface)
                if not edge_available:
                    checked_ifs.append((chain1, chain2, source))
                    break
                elif not xl_ok:
                    break
            elif chain2 not in self.groups:
                edge_available, xl_ok = self.clash_append(chain2, chain1, interface)
                if not edge_available:
                    checked_ifs.append((chain1, chain2, source))
                    break
                elif not xl_ok:
                    break
            elif chain1 not in self.groups[chain2]:
                edge_available, xl_ok = self.clash_merge(chain1, chain2, interface)
                if not edge_available:
                    checked_ifs.append((chain1, chain2, source))
                    break
                elif not xl_ok:
                    break
        return checked_ifs


# define a sampling process
def assembly_process_xl(interfaces, assembler: XLink_Assembler, max_retries: int = 5):
    assembler.reset()
    stime = datetime.now()
    assemblied_all, xl_excluded_ifs = assembler.assembly_in_order(interfaces)
    n_retry = 0
    while not assemblied_all and xl_excluded_ifs and n_retry < max_retries:
        prev_xl_excluded_ifs = xl_excluded_ifs.copy()
        new_ifs = utils.create_path_df(interfaces, xl_excluded_ifs)
        assemblied_all, xl_excluded_ifs = assembler.assembly_in_order(new_ifs)
        xl_same = (len(xl_excluded_ifs) == len(prev_xl_excluded_ifs) and
                   all(a.equals(b) for a, b in zip(xl_excluded_ifs, prev_xl_excluded_ifs)))
        if xl_same:
            break
        n_retry += 1

    if assemblied_all:
        assembler.orientation = {
            chain_id: assembler.orientation[chain_id] for chain_id in config.chain_ids
        }
        chain_id = config.chain_ids[0]
        if_score = assembler.calculate_interface_score(chain_id)
        xl_satisfaction = utils.calculate_xl_satisfaction(config.xl_total_data, assembler.xl_ok_data)
        new_if_score = xl_satisfaction * if_score
        assemblied_result = {
            "if_score": new_if_score, "path": assembler.path, "orientation": assembler.orientation,
            "xl_satisfaction": xl_satisfaction, "original_if_score": if_score
        }
        etime = datetime.now()
        run_time = etime - stime
        return True, True, assemblied_result, run_time
    else:
        length = 0
        chain = ""
        group_length = 0
        for chain_id in config.chain_ids:
            if chain_id in assembler.groups:
                group_length += len(assembler.groups[chain_id])
                if len(assembler.groups[chain_id]) > length:
                    length = len(assembler.groups[chain_id])
                    chain = chain_id
            else:
                group_length += 1
        mean_length = group_length / config.chain_num
        if_score = assembler.calculate_interface_score(chain)
        total_xls, satisfied_xls = assembler.count_group_xls(assembler.groups[chain], {}, {})
        xl_satisfaction = utils.calculate_xl_satisfaction(total_xls, satisfied_xls)
        new_if_score = if_score * xl_satisfaction
        subunits_connectivity = True
        checked_chain_ids = set()
        for chain_id in config.chain_ids:
            if chain_id in checked_chain_ids:
                continue
            if chain_id in config.subunits_id_mapping:
                full_chain_id = config.subunits_id_mapping[chain_id]
                subunits = set(config.subunits_order[full_chain_id].keys())
                group_set = set(assembler.groups.get(chain_id, []))
                if not subunits.issubset(group_set):
                    subunits_connectivity = False
                    break
                checked_chain_ids.update(subunits)
        assemblied_result = {
            "if_score": new_if_score, "max_length": length, "mean_length": mean_length,
            "xl_satisfaction": xl_satisfaction, "original_if_score": if_score,
            "groups": assembler.groups, "orientation": assembler.orientation
        }
        etime = datetime.now()
        run_time = etime - stime
        return False, subunits_connectivity, assemblied_result, run_time

# define a sampling task
def assembly_task_xl(interfaces, xl_threshold: float = 0.7):
    try:
        single_assembler = XLink_Assembler(xl_threshold)
        assemblied_all, subunits_connectivity, assemblied_result, run_time = assembly_process_xl(
            interfaces, single_assembler, 5
        )
        del single_assembler
        gc.collect()
        return assemblied_all, subunits_connectivity, assemblied_result, run_time
    except Exception as e:
        print(f"[ERROR] Error in assembly sampling process: {e}", flush=True)
        gc.collect()
        return None, None, None, None

# define the initial sampling stage
def generate_initial_population_xl(num_cpus, total_tasks, fout, xl_threshold: float = 0.7):
    successful_results = []
    failed_results_connected = []
    failed_results = []
    stime = datetime.now()
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = []
        for i in range(total_tasks):
            if i == 0:
                interfaces = utils.sort_dataframe(
                    config.interfaces, "score", "scaled_score2", "chain1", "chain2"
                )
            elif i == 1:
                interfaces = utils.sort_dataframe(
                    config.interfaces, "mean_score", "score", "chain1", "chain2"
                )
            elif i == 2:
                interfaces = utils.sort_dataframe(
                    config.interfaces, "scaled_score2", "score", "chain1", "chain2"
                )
            else:
                interfaces = utils.shuffle_dataframe(config.interfaces, config.RANDOM_SEED + i)
            futures.append(executor.submit(assembly_task_xl, interfaces, xl_threshold))
        num = 0
        for future in as_completed(futures):
            num += 1
            assemblied_all, subunits_connectivity, assemblied_result, run_time = future.result()
            print(f"Round {num} in initial random sampling cost {run_time}",
                  file=fout, flush=True)
            if assemblied_all:
                successful_results.append(assemblied_result)
            elif subunits_connectivity:
                failed_results_connected.append(assemblied_result)
            elif assemblied_result:
                failed_results.append(assemblied_result)
    gc.collect()
    etime = datetime.now()
    print(f"\nGenerate initial population cost {etime - stime}", file=fout, flush=True)
    print(f"\nGenerate initial population cost {etime - stime}", flush=True)
    if successful_results:
        return True, True, successful_results
    elif failed_results_connected:
        return False, True, failed_results_connected
    else:
        return False, False, failed_results

