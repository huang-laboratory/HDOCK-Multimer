import re
import os
import shutil
import subprocess
import numpy as np
from chainlist import ChainList
from multi_utils import (
    DEFAULT_MODEL_FOLDER,
    PDB_PARSER,
    get_bb_atoms,
    get_model_file,
    neighbor_atoms_number,
    get_complex_name,
    inter_atoms,
    get_reverse_model
)
from typing import Dict, List
import sys 


SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
# sys.path.append(BASE_PATH)
transchain = os.path.join(TOOL_PATH, "transchain")


class Picker:
    """ Factory class for picker """
    chains: ChainList
    scores: np.ndarray
    chains_poses: List[str]
    complexes: List[List[int]]
    merged_chains: Dict[int, int]

    def __init__(self, fout, normalized: bool, model: int, append: bool = True, model_path: str = DEFAULT_MODEL_FOLDER):
        self.fout = fout
        self.normalized = normalized
        self.model = model
        self.append = append
        self.model_path = model_path

    def __len__(self):
        return len(self.chains)

    def _get_scores(self):
        """ parse pairwise docking scores, save to self.scores and print """
        chains = self.chains
        chain_num = len(chains)
        self.scores = np.zeros((chain_num, chain_num))
        for i, _ in enumerate(chains):
            for j in range(i):
                self.scores[i][j] = self.scores[j][i]
            for j in range(i+1, chain_num):
                model_filename = self.chains.get_model_file(i, j, 1, self.model_path)
                with open(model_filename, 'r', encoding='utf-8') as model_file:
                    for line in model_file:
                        if "REMARK Score" in line:
                            self.scores[i][j] = float(
                                re.match(r"REMARK Score:\s*(-[.0-9]+)", line).group(1)
                            )
                            break
                    if self.normalized:
                        lig_atoms = list(
                            PDB_PARSER.get_structure("lig", model_file).get_atoms()
                        )
                        rec_atoms = list(
                            PDB_PARSER.get_structure("rec", chains[i]).get_atoms()
                        )
                        inter_rec, inter_lig = inter_atoms(rec_atoms, lig_atoms)
                        # self.scores[i][j] /= len(inter_rec) * len(inter_lig)
                        self.scores[i][j] /= len(rec_atoms) * len(lig_atoms)
        np.set_printoptions(precision=2, linewidth=10 * (chain_num + 1), suppress=False)
        print(self.scores, file=self.fout, flush=True)

    def _new_complex(self, chain_idx: int, partner_idx: int):
        """
        pick a new complex
        append a list of interacting index to self.complexes
        update self.merged_chains and corresponding chains_poses

        Args:
            - chain_idx (int): receptor index in self.chains, fixed in docking
            - partner_idx (int): ligand index in self.chains, moving in docking
        """
        comp_idx = len(self.complexes)  # new complex imdex
        self.merged_chains.update([(chain_idx, comp_idx), (partner_idx, comp_idx)])
        self.complexes.append([chain_idx, partner_idx])  # update complex information
        self.chains_poses[chain_idx] = self.chains[chain_idx]  # fix receptor
        pose = get_model_file(
            self.chains[chain_idx], self.chains[partner_idx], 1, self.model_path
        )  # transformed ligand file
        init_model = self.chains.get_model_file(
            chain_idx, partner_idx, 1, self.model_path
        )  # initial transformed ligand file, corresponding to repeat situations
        if pose != init_model:
            shutil.copy(init_model, pose)
        self.chains_poses[partner_idx] = pose
        self.chains.append_trans(
            partner_idx, self.chains.get_outfile_name(chain_idx, partner_idx), False
        )
        score = self.scores[chain_idx][partner_idx]
        print(
            f"[INFO] complex {comp_idx}, {self.chains[chain_idx]} {self.chains[partner_idx]} with {score}",
            file=self.fout, flush=True
        )

    def _append(self, rec_idx: int, lig_idx: int, high=False):
        """
        append a chain to another bound chain

        Args:
            rec_idx (int): bound chain index
            lig_idx (int): unbound chain index, which is to be appended
            high (bool): pairwise score is too high, default is False
        """
        comp_idx = self.merged_chains[rec_idx]
        receptor = self.chains[rec_idx]
        ligand = self.chains[lig_idx]

        # make sure the moving chain is smaller
        if rec_idx > lig_idx:
            model_file = self.chains.get_reverse_model(
                self.chains.init_chain_filename(receptor), self.chains.init_chain_filename(ligand), 1, self.model_path
            )
        else:
            model_file = get_model_file(receptor, ligand, 1, self.model_path)
            init_model = self.chains.get_model_file(receptor, ligand, 1, self.model_path)
            if model_file != init_model:
                shutil.copy(init_model, model_file)

        # apply receptor's transformation to the ligand model
        # the model file will be used only once, so transform it inplace
        for trans, inverse in self.chains.transforms[rec_idx]:
            if inverse:
                subprocess.run(f"{transchain} {model_file} {trans} -r", shell=True, check=True)
            else:
                subprocess.run(f"{transchain} {model_file} {trans}", shell=True, check=True)

        # get backbone atoms of receptor, ligand, and the complex excluding the receptor
        rec_atoms = get_bb_atoms(PDB_PARSER.get_structure("chain", receptor))
        lig_atoms = get_bb_atoms(PDB_PARSER.get_structure("chain", model_file))
        other_atoms = []
        for i in self.complexes[comp_idx]:
            if i != rec_idx:
                other_atoms += list(
                    get_bb_atoms(PDB_PARSER.get_structure("chain", self.chains_poses[i]))
                )
        # append successfully if additional contacts in tolerance
        additional_atoms = neighbor_atoms_number(other_atoms, lig_atoms)
        # inter_rec, inter_lig = inter_atoms(rec_atoms, lig_atoms)
        # tolerance = (len(inter_rec) + len(inter_lig)) * 0.1
        tolerance = 30
        if additional_atoms < tolerance:
            score = self.scores[rec_idx][lig_idx]
            if high:
                print(f"[WARN] append {ligand} to {comp_idx} bust {score} too high", file=self.fout)
            else:
                self.merged_chains[lig_idx] = comp_idx
                self.complexes[comp_idx].append(lig_idx)
                self.chains_poses[lig_idx] = model_file
                for trans, inverse in self.chains.transforms[rec_idx]:
                    self.chains.append_trans(lig_idx, trans, inverse)
                if rec_idx > lig_idx:
                    self.chains.append_trans(lig_idx, self.chains.get_outfile_name(lig_idx, rec_idx), True)
                else:
                    self.chains.append_trans(lig_idx, self.chains.get_outfile_name(rec_idx, lig_idx), False)
            print(
                f"[INFO] append {ligand} to {comp_idx} contact to {receptor} with {score}",
                file=self.fout, flush=True
            )
        else:
            print(
                f"[WARN] {ligand} additional {additional_atoms} contact interactions with {receptor}",
                file=self.fout
            )

    def _merge(self):
        """ merge each picked chains group and collect unpicked chains """
        results: List[str] = []
        for comp in self.complexes:
            filenames = [self.chains[i] for i in comp]
            merged_filename = self.model_path + get_complex_name(*filenames) + ".pdb"
            pose_filenames = [self.chains_poses[i] for i in comp]
            with open(merged_filename, 'w+', encoding='utf-8') as fmerged:
                for pose_filename in pose_filenames:
                    chain_id: str = " "
                    new_id: str = " "
                    with open(pose_filename, 'r', encoding='utf-8') as fpose:
                        for line in fpose:
                            if line.startswith("TER"):
                                fmerged.write(line)
                            if line.startswith("ATOM"):
                                if chain_id != line[21]:
                                    chain_id = line[21]
                                    if chain_id in self.chains.chain_ids:
                                        new_id = self.chains.unique_chain_id()
                                    else:
                                        new_id = chain_id
                                        self.chains.chain_ids.add(chain_id)
                                fmerged.write(line[:21] + new_id + line[22:])
            results.append(merged_filename)

        for i, pose in enumerate(self.chains_poses):
            if pose == "":
                results.append(self.chains[i])

        return ChainList(results)

    def pick(self, chainlist: ChainList):
        """ super pick method for all pickers """
        self.chains = chainlist
        self.chains_poses = ["" for _ in self.chains]
        self.complexes: List[List[int]] = []
        self.merged_chains: Dict[int, int] = {}

    @staticmethod
    def get_picker(fout, normalized: bool, method: int, model: int, append: bool, model_path: str = DEFAULT_MODEL_FOLDER):
        """
        factory method for picker
        Args:
            normalized (bool): whether to normalize score by interacting atoms
            model_path (str): path of models generated by docking
            method (int): pick method id
        Raises:
            ValueError: method id is not in one of [0, 1, 2, 3]
        Returns:
            Picker subclass instance: _description_
        """
        if method == 0:
            return _PickerEvery(fout, normalized, model, append, model_path)
        if method == 1:
            return _PickerOnly(fout, normalized, model, append, model_path)
        if method == 2:
            return _BasePicker(fout, normalized, model, append, model_path)
        if method == 3:
            return _HigherPicker(fout, normalized, model, append, model_path)
        raise ValueError(f"[ERROR] Invalid method id: {method}, which should be one of [0,1,2,3]")


class _PickerEvery(Picker):
    """ pick every chain partner with the highest score by the given order """
    def pick(self, chainlist: ChainList):
        super().pick(chainlist)
        self._get_scores()
        min_idx_list = self.scores.argmin(1)
        print(min_idx_list, file=self.fout)

        for chain_idx, partner_idx in enumerate(min_idx_list):
            if chain_idx not in self.merged_chains:
                if partner_idx not in self.merged_chains:
                    # new complex
                    if chain_idx > partner_idx:
                        chain_idx, partner_idx = partner_idx, chain_idx
                    self._new_complex(chain_idx, partner_idx)
                else:
                    # partner is not fixed which may cause additional interactions
                    if self.chains_poses[partner_idx] != chainlist[partner_idx]:
                        continue
                    # add chain to previous complex if there is no clash
                    self._append(partner_idx, chain_idx)

        return self._merge()


class _PickerOnly(Picker):
    """ only pick the complex with the highest score in each iteration """
    def pick(self, chainlist: ChainList):
        super().pick(chainlist)
        self._get_scores()
        chain_idx, partner_idx = np.unravel_index(self.scores.argmin(), self.scores.shape)
        self._new_complex(int(chain_idx), int(partner_idx))
        return self._merge()


class _BasePicker(Picker):
    """ sort the lowest energies for each chain in ascending, pick the base chain by order """
    def pick(self, chainlist: ChainList, cutoff: float = None):
        super().pick(chainlist)
        self._get_scores()
        # list of partner chain indices with the minimal score
        min_score_list = self.scores.min(1)
        min_idx_list = [[] for _ in min_score_list]
        for i, min_value in enumerate(min_score_list):
            for j, value in enumerate(self.scores[i]):
                if value == min_value:
                    min_idx_list[i].append(j)
        min_list = np.sort(
            np.array(
                [
                    (chain_idx, partner_idx, self.scores[chain_idx][partner_idx])
                    for chain_idx, partner_indices in enumerate(min_idx_list)
                    for partner_idx in partner_indices
                ],
                dtype=[("chain_idx", int), ("partner_idx", int), ("min_score", float)]
            ), kind="stable", order="min_score"
        )
        for chain_idx, partner_idx, min_score in min_list:
            if (
                cutoff is not None
                and (min_list["min_score"][0] - min_score) / min_list["min_score"][0] > cutoff
            ):
                break

            if chain_idx not in self.merged_chains:
                if partner_idx not in self.merged_chains:
                    if chain_idx > partner_idx:
                        print("[WARN] chain index larger", file=self.fout)
                        continue
                    self._new_complex(chain_idx, partner_idx)
                elif self.append:
                    # partner is not fixed which may cause additional interactions
                    # if self.chains_poses[partner_index] != chainlist[partner_index]:
                    #     print("[WARN] partner is not fixed, file=self.fout")
                    #     continue

                    # add chain to previous complex if there is no clash and score is not too high
                    comp_idx = self.complexes[self.merged_chains[partner_idx]]
                    comp_min_score = self.scores[comp_idx[0], comp_idx[1]]
                    high = min_score > (comp_min_score // -200) * -100
                    self._append(partner_idx, chain_idx, high)
        return self._merge()


class _HigherPicker(Picker):
    """ sort every energies and pick the higher ones to generate comple """
    def pick(self, chainlist: ChainList, cutoff: float = None):
        super().pick(chainlist)
        self._get_scores()
        min_list = np.sort(
            np.array(
                [(x, y, element) for (x, y), element in np.ndenumerate(self.scores)],
                dtype=[("chain_idx", int), ("partner_idx", int), ("min_score", float)]
            ), order="min_score"
        )
        for chain_idx, partner_idx, min_score in min_list:
            if chain_idx == partner_idx:
                continue
            if (
                cutoff is not None
                and (min_list["min_score"][0] - min_score) / min_list["min_score"][0] > cutoff
            ):
                break

            if chain_idx not in self.merged_chains:
                if partner_idx not in self.merged_chains:
                    if chain_idx > partner_idx:
                        print("[WARN] chain index larger", file=self.fout)
                        continue
                    self._new_complex(chain_idx, partner_idx)
                elif self.append:
                    # partner is not fixed which may cause additional interactions
                    # if self.chains_poses[partner_index] != chainlist[partner_index]:
                    #     print("[WARN] partner is not fixed", file=self.fout)
                    #     continue

                    # add chain to previous complex if there is no clash
                    self._append(partner_idx, chain_idx)
        return self._merge()
