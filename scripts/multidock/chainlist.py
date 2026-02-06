"""
ChainList implementation module,
which maintains the list of chains fro each docking round
"""
import shutil
import os
import subprocess
import sys 
from typing import Union, Iterable, List, Dict
import multi_utils
from transform import Transform
from multi_utils import get_atoms_num

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
sys.path.append(BASE_PATH)
transchain = os.path.join(TOOL_PATH, "transchain")


class ChainList:
    """ List of chains wrapper """
    def __init__(self, init_chains: Iterable[str], repeat: Iterable[int] = None) -> None:
        # valid if all chain files containing ATOM
        for i in init_chains:
            empty = True
            with open(i, 'r', encoding='utf-8') as fchain:
                for line in fchain:
                    if line.startswith("ATOM"):
                        empty = False
                        break
                if empty:
                    raise OSError(f"Empty file: {i}")
        if repeat is None:
            repeat = [1 for _ in init_chains]

        self.num = sum(repeat)  # number of all chains
        self.init_chains, self.repeat = zip(
            *sorted(
                zip(init_chains, repeat),
                key=lambda pair: get_atoms_num(pair[0]),
                reverse=True
            )
        )  # files with more atoms are placed at the top
        self.unique: Dict[str, int] = {}
        self.chains: List[str] = []  # save filenames for all chains
        self.chain_ids = set()
        self.transforms: List[Transform] = []

    def __len__(self):
        return self.num

    def __iter__(self):
        return iter(self.chains)

    def __getitem__(self, key) -> str:
        return self.chains[key]

    def __str__(self) -> str:
        return str(self.chains)

    def duplicate(self):
        """ duplicate repeated chains, should be called after changing work dir """
        for i, chain in enumerate(self.init_chains):
            if self.repeat[i] == 1:
                self.unique[chain] = i
                self.chains.append(chain)
                continue

            for j in range(self.repeat[i]):
                duplication = chain[:-4] + str(j + 1) + ".pdb"
                shutil.copy(chain, duplication)
                self.unique[chain] = i
                self.unique[duplication] = i
                self.chains.append(duplication)

        self.transforms = [Transform() for _ in self.chains]

    def chain_repeat(self, chain: Union[str, int]):
        """ get repeat number of chain """
        if isinstance(chain, str):
            return self.repeat[self.unique[chain]]
        return self.repeat[self.unique[self.chains[chain]]]

    def init_chain_filename(self, chain: Union[str, int]):
        """ get initial chain filename """
        if isinstance(chain, str):
            return self.init_chains[self.unique[chain]]
        return self.init_chains[self.unique[self.chains[chain]]]

    def get_model_file(self, rec_idx: Union[str, int], lig_idx: Union[str, int],
                       model_number: int, folder: str):
        """ get model file name """
        return multi_utils.get_model_file(
            self.init_chain_filename(rec_idx), self.init_chain_filename(lig_idx), model_number, folder
        )

    def get_reverse_model(self, receptor: str, ligand: str,
                          model_number: int, folder=multi_utils.DEFAULT_MODEL_FOLDER):
        """ reverse positions of receptor and ligand """
        output_file = multi_utils.get_outfile_name(
            self.init_chain_filename(ligand), self.init_chain_filename(receptor)
        )
        reverse_model = multi_utils.get_model_file(
            receptor, ligand, model_number, folder
        )
        shutil.copy(ligand, reverse_model)
        subprocess.run(f"{transchain} {reverse_model} {output_file} -r", shell=True, check=True)
        return reverse_model

    def get_outfile_name(self, rec_idx: Union[str, int], lig_idx: Union[str, int]):
        """ get hdock output filename """
        return multi_utils.get_outfile_name(
            self.init_chain_filename(rec_idx), self.init_chain_filename(lig_idx)
        )

    def unique_chain_id(self):
        """ generate new chain ID """
        new_ids = ([chr(i) for i in range(ord("A"), ord("Z") + 1)] +
                   [chr(i) for i in range(ord("a"), ord("z") + 1)] +
                   [chr(i) for i in range(ord("0"), ord("9") + 1)])
        for new_id in new_ids:
            if new_id not in self.chain_ids:
                self.chain_ids.add(new_id)
                return new_id
        return "+"

    def append_trans(self, idx: int, file: str, inverse: bool):
        """ record new docking transformation """
        self.transforms[idx].append(file, inverse)
