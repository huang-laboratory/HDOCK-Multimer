""" Some configs and litter functions for multi docking """
from itertools import combinations
import re
import os
import shutil
import subprocess
from pathlib import Path
import argparse
from typing import Generator, Iterable, Union, List, Dict
import sys 

SCRIPT_PATH = os.path.dirname(os.path.abspath(__file__))
BASE_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
TOOL_PATH = os.path.abspath(os.path.join(BASE_PATH, "..", "tools"))
sys.path.append(BASE_PATH)
transchain = os.path.join(TOOL_PATH, "transchain")

from parse_pdb import Atom, Residue, Chain, Structure
from parse_pdb import PDBParser, KDNeighborSearch


DEFAULT_MODEL_FOLDER = "models/"
DEFAULT_OUTFILE_PATH = "pairwise_out/"
PDB_PARSER = PDBParser()


def get_atoms_num(file: str):
    """ Return number of atoms from the file """
    with open(file, 'r', encoding='utf-8') as fpdb:
        return sum(1 for line in fpdb if line.startswith("ATOM "))

def get_ca_atoms(entity: Union[Structure, Chain, Residue]) -> Generator[Atom, None, None]:
    """ Return CA atoms """
    return (atom for atom in entity.get_atoms() if atom.name == "CA")

def get_bb_atoms(entity: Union[Structure, Chain, Residue]):
    """ Return backbone atoms """
    return (atom for atom in entity.get_atoms() if atom.name in ["CA", "C", "N", "O"])

def get_complex_name(*args: str) -> str:
    """ Trim chain file path and concatenate to a complex name """
    patten = r"(.*[/.]|^)(.+).pdb"
    return "".join([re.match(patten, i).group(2) for i in args])

def get_outfile_name(receptor: str, ligand: str):
    """ Return HDOCK output filename, like reclig.out """
    return get_complex_name(receptor, ligand) + ".out"

def get_model_file(receptor: str, ligand: str, model_number: int, folder=DEFAULT_MODEL_FOLDER):
    """ Return the model filename, like folder/rec.lig.pdb """
    # return f"{folder}{get_complex_name(receptor)}.{model_number}{get_complex_name(ligand)}.pdb"
    return f"{folder}{get_complex_name(receptor)}.{get_complex_name(ligand)}.pdb"

def get_reverse_model(receptor: str, ligand: str, model_number: int, folder=DEFAULT_MODEL_FOLDER):
    """ Fix ligand and move receptor reverse """
    output_file = get_outfile_name(ligand, receptor)
    reverse_model = get_model_file(receptor, ligand, model_number, folder)
    shutil.copy(ligand, reverse_model)
    subprocess.run(f"{transchain} {reverse_model} {output_file} -r", shell=True, check=True)
    return reverse_model

def get_model_file_ex(chain1: str, chain2: str, model_number: int, folder=DEFAULT_MODEL_FOLDER):
    """ Check whether the model file exists """
    if Path(model := get_model_file(chain1, chain2, model_number, folder)).exists():
        return chain1, chain2, model
    if Path(model := get_model_file(chain2, chain1, model_number, folder)).exists():
        return chain2, chain1, model
    raise FileNotFoundError(f"model of {chain1} and {chain2} is not found.")

def neighbor_atoms(query_atoms: List[Atom], target_atoms: Iterable[Atom], radius: float = 10.0):
    ns = KDNeighborSearch(query_atoms)
    return (neighbor for atom in target_atoms for neighbor in ns.search(atom.coord, radius))

def neighbor_atoms_number(query_atoms: List[Atom], target_atoms: Iterable[Atom], radius: float = 10.0):
    return len(list(neighbor_atoms(query_atoms, target_atoms, radius)))

def atoms_are_contacted(query_atoms: List[Atom], target_atoms: Iterable[Atom], radius: float = 10.0):
    return next(neighbor_atoms(query_atoms, target_atoms, radius), None) is not None

def entities_are_contacted(
    entity1: Union[Structure, Chain, Residue],
    entity2: Union[Structure, Chain, Residue],
    radius: float = 10.0
):
    return (
        next(
            neighbor_atoms(list(entity1.get_atoms()), entity2.get_atoms(), radius), None
        ) is not None
    )

def inter_res(rec_atoms: List[Atom], lig_atoms: Iterable[Atom], radius: float = 10.0):
    ns = KDNeighborSearch(rec_atoms)
    inter_rec, inter_lig = set(), set()
    for atom in lig_atoms:
        neighbors = ns.search(atom.coord, radius)
        if neighbors:
            inter_lig.add(atom.get_parent())
            inter_rec.update(neighbor.get_parent() for neighbor in neighbors)
    return inter_rec, inter_lig

def inter_atoms(rec_atoms: List[Atom], lig_atoms: Iterable[Atom], radius: float = 10.0):
    ns = KDNeighborSearch(rec_atoms)
    inter_rec, inter_lig = set(), set()
    for atom in lig_atoms:
        neighbors = ns.search(atom.coord, radius)
        if neighbors:
            inter_lig.add(atom)
            inter_rec.update(neighbors)
    return inter_rec, inter_lig

def get_interfaces(chains: Iterable[Chain], radius: float = 10.0):
    return (
        (chain1, chain2) for chain1, chain2 in combinations(chains, 2)
        if entities_are_contacted(chain1, chain2, radius)
    )

def merge(chains_filenames: List[str], outfolder: str=DEFAULT_MODEL_FOLDER):
    """
    Merge all chains to a new pdb file by the list order
    Arguments:
        - chains_filenames: list of filenames of chains
        - outfolder: path to the new pdb file, default is current directory

    Returns:
        Path to the new pdb file. If chains_filenames length is 1, return the file
    """
    if len(chains_filenames) == 1:
        return chains_filenames[0]
    merged_filename = outfolder + get_complex_name(*chains_filenames) + ".pdb"
    with open(merged_filename, 'w+', buffering=1024 * 10, encoding='utf-8') as fmerge:
        for chain_filename in chains_filenames:
            with open(chain_filename, 'r', encoding='utf-8') as fchain:
                for line in fchain:
                    if line.startswith(("ATOM", "TER")):
                        fmerge.write(line)
    return merged_filename

def pdb_file_checker(filename: str):
    """ Checker for argparse, check if filename is a valid PDB file """
    if not os.path.isfile(filename):
        raise argparse.ArgumentTypeError('"' + filename + '" is not a file')
    if not filename.endswith("pdb"):
        raise argparse.ArgumentTypeError("argument filename must be of type *.pdb")
    return filename
