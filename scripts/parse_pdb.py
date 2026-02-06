import numpy as np
import copy
from scipy.spatial import KDTree


class Atom:
    def __init__(self, atm_name, alt_loc, coord, occupancy, bfactor, element):
        self.name = atm_name
        self.alt_loc = alt_loc
        self.coord = np.array(coord)
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.element = element
        self._parent = None

    def get_name(self):
        return self.name

    def get_coord(self):
        return self.coord

    def get_bfactor(self):
        return self.bfactor

    def transform(self, rot, trans):
        self.coord = np.dot(self.coord, rot) + trans

    def get_parent(self):
        return self._parent


class Residue:
    def __init__(self, res_name, resseq, icode):
        self.resname = res_name
        self.resseq = resseq
        self.icode = icode
        self.atoms = {}
        self._parent = None

    def __getitem__(self, atm_name):
        if atm_name in self.atoms:
            return self.atoms[atm_name]
        raise KeyError(f"Atom {atm_name} not found!")

    def add_atom(self, atom):
        atom._parent = self
        self.atoms[atom.name] = atom

    def get_atoms(self):
        return self.atoms.values()

    def get_parent(self):
        return self._parent


class Chain:
    def __init__(self, chain_id):
        self.id = chain_id
        self.residues = []
        self._parent = None

    def add_residue(self, residue):
        residue._parent = self
        self.residues.append(residue)

    def get_residues(self):
        return self.residues

    def get_atoms(self):
        for residue in self.residues:
            yield from residue.get_atoms()
 
    def copy(self):
        return copy.deepcopy(self)

    def get_parent(self):
        return self._parent


class Model:
    def __init__(self, model_id):
        self.id = model_id
        self.chains = {}
        self._parent = None

    def __getitem__(self, chain_id):
        if chain_id in self.chains:
            return self.chains[chain_id]
        raise KeyError(f"Chain {chain_id} not found!")

    def add_chain(self, chain):
        chain._parent = self
        self.chains[chain.id] = chain

    def get_chains(self):
        return self.chains.values()

    def get_atoms(self):
        for chain in self.get_chains():
            yield from chain.get_atoms()

    def get_parent(self):
        return self._parent


class Structure:
    def __init__(self, struct_id):
        self.id = struct_id
        self.models = []

    def __getitem__(self, model_id):
        if len(self.models) <= model_id:
            raise IndexError(f"Model index {model_id} out of range!")
        else:
            return self.models[model_id]

    def add_model(self, model):
        model._parent = self
        self.models.append(model)

    def get_models(self):
        return self.models

    def get_atoms(self):
        for model in self.models:
            yield from model.get_atoms()

    def get_chains(self):
        for model in self.models:
            yield from model.get_chains()


class PDBParser:
    def __init__(self):
        self.structure = None
        self.model = None
        self.chain = None
        self.residue = None

    def get_structure(self, struct_id, pdb_file):
        self.reset()
        self.structure = Structure(struct_id)
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith("MODEL"):
                    model_id = int(line[10:14].strip())
                    self.model = Model(model_id)
                    self.structure.add_model(self.model)
                elif line.startswith("ENDMDL"):
                    self.model = None
                #elif line.startswith("ATOM") or line.startswith("HETATM"):
                elif line.startswith("ATOM"):
                    #self._parse_atom_line(line)
                    atom_name = line[12:16].strip()
                    if not atom_name.startswith("H"):
                        self._parse_atom_line(line)
                elif line.startswith("END"):
                    break
        return self.structure

    def _parse_atom_line(self, line):
        atm_name = line[12:16].strip()
        alt_loc = line[16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21].strip()
        res_seq = int(line[22:26].strip())
        icode = line[26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        if line[54:60].strip() == "":
            occupancy = ""
            bfactor = ""
            element = ""
        else:
            occupancy = float(line[54:60].strip())
            bfactor = float(line[60:66].strip())
            element = line[76:78].strip()

        atom = Atom(atm_name, alt_loc, (x, y, z), occupancy, bfactor, element)
        if self.model is None:
            self.model = Model(0)
            self.structure.add_model(self.model)
        if (self.chain is None) or (self.chain.id != chain_id):
            self.chain = Chain(chain_id)
            self.model.add_chain(self.chain)
        if (self.residue is None) or (self.residue.resseq != res_seq) or (self.residue.icode != icode):
            self.residue = Residue(res_name, res_seq, icode)
            self.chain.add_residue(self.residue)
        self.residue.add_atom(atom)

    def reset(self):
        self.structure = None
        self.model = None
        self.chain = None
        self.residue = None


class KDNeighborSearch:
    def __init__(self, atoms):
        self.atoms = atoms
        self.coords = np.array([atom.coord for atom in atoms])
        self.kdtree = KDTree(self.coords)

    def search(self, center_coord, radius):
        indices = self.kdtree.query_ball_point(center_coord, radius)
        return [self.atoms[i] for i in indices]

    def search_atom(self, atom, radius):
        return self.search(atom.coord, radius)

    def search_all_pairs(self, radius):
        pairs = self.kdtree.query_pairs(radius)
        return [(self.atoms[i], self.atoms[j]) for i, j in pairs]


def load_format_line(format_file):
    format_lines = {}
    with open(format_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atm_name = line[12:16].strip()
                res_name = line[17:20].strip()
                format_lines.setdefault(res_name, {})[atm_name] = line
    return format_lines


def format_atom_line(natm, line, atom: Atom, residue: Residue, chain: Chain) -> str:
    if atom.get_bfactor() == "":
        line = (
            line[:6] + f"{natm:5d}" + line[11:16] + f"{atom.alt_loc:1}" +
            line[17:21] + chain.id + f"{residue.resseq:4d}" + f"{residue.icode:1}" +
            line[27:30] + f"{atom.coord[0]:8.3f}" + f"{atom.coord[1]:8.3f}" + f"{atom.coord[2]:8.3f}" +
            f"\n"
        )
    else:
        line = (
            line[:6] + f"{natm:5d}" + line[11:16] + f"{atom.alt_loc:1}" +
            line[17:21] + chain.id + f"{residue.resseq:4d}" + f"{residue.icode:1}" +
            line[27:30] + f"{atom.coord[0]:8.3f}" + f"{atom.coord[1]:8.3f}" + f"{atom.coord[2]:8.3f}" +
            f"{atom.occupancy:6.2f}" + f"{atom.bfactor:6.2f}" + line[66:]
        )
    return line


def write_pdb_file(data, pdb_file, format_lines):
    with open(pdb_file, 'w') as file:
        if isinstance(data, Structure):
            for model in data.get_models():
                natm = 0
                for chain in model.get_chains():
                    for residue in chain.get_residues():
                        if residue.resname not in format_lines:
                            continue
                        for atom in residue.get_atoms():
                            natm += 1
                            format_line = format_lines[residue.resname][atom.name]
                            file.write(format_atom_line(natm, format_line, atom, residue, chain))
                    file.write("TER\n")
        elif isinstance(data, Model):
            natm = 0
            for chain in data.get_chains():
                for residue in chain.get_residues():
                    if residue.resname not in format_lines:
                        continue
                    for atom in residue.get_atoms():
                        natm += 1
                        format_line = format_lines[residue.resname][atom.name]
                        file.write(format_atom_line(natm, format_line, atom, residue, chain))
                file.write("TER\n")
        elif isinstance(data, Chain):
            natm = 0
            for residue in data.get_residues():
                if residue.resname not in format_lines:
                    continue
                for atom in residue.get_atoms():
                    natm += 1
                    format_line = format_lines[residue.resname][atom.name]
                    file.write(format_atom_line(natm, format_line, atom, residue, data))
            file.write("TER\n")
        elif isinstance(data, dict):
            natm = 0
            for chain in data.values():
                for residue in chain.get_residues():
                    if residue.resname not in format_lines:
                        continue
                    for atom in residue.get_atoms():
                        natm += 1
                        format_line = format_lines[residue.resname][atom.name]
                        file.write(format_atom_line(natm, format_line, atom, residue, chain))
                file.write("TER\n")
        elif isinstance(data, list):
            natm = 0
            for chain in data:
                for residue in chain.get_residues():
                    if residue.resname not in format_lines:
                        continue
                    for atom in residue.get_atoms():
                        natm += 1
                        format_line = format_lines[residue.resname][atom.name]
                        file.write(format_atom_line(natm, format_line, atom, residue, chain))
                file.write("TER\n")
