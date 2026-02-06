import sys
from openmm import unit, Platform, CustomExternalForce, NonbondedForce, VerletIntegrator
from openmm import LocalEnergyMinimizer
from openmm.app import ForceField, Simulation, HBonds, NoCutoff, Modeller
from openmm.app import PDBFile
from openmm.app import element as elem
from pdbfixer import PDBFixer


STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEP", "TPO", "PTR"
}


def kcalA2_to_kJnm2(k):
    # 1 kcal/mol/Å^2 = 418.4 kJ/mol/nm^2
    return k * 418.4

# load force field parameters
def try_forcefield_with_phos():
    ff_xmls = ["amber14/protein.ff14SB.xml"]
    for cand in ["amber/phosaa10.xml", "amber14/phosaa10.xml", "amber10/phosaa10.xml"]:
        try:
            ForceField(*ff_xmls, cand)
            ff_xmls.append(cand)
            break
        except Exception:
            pass
    return ForceField(*ff_xmls)

# add position constraints on backbone atoms (CA, C, N)
def add_backbone_position_restraints(system, topology, positions, k_kcal_per_A2=1.0):
    k = kcalA2_to_kJnm2(k_kcal_per_A2)
    restr = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restr.addGlobalParameter("k", k)
    restr.addPerParticleParameter("x0")
    restr.addPerParticleParameter("y0")
    restr.addPerParticleParameter("z0")

    for atom in topology.atoms():
        if atom.residue.name not in STANDARD_AA:
            continue
        if atom.name.strip() in ("CA", "C", "N"):
            xyz_nm = positions[atom.index].value_in_unit(unit.nanometer)
            restr.addParticle(atom.index, [xyz_nm[0], xyz_nm[1], xyz_nm[2]])
    system.addForce(restr)

# remove the hydrogen atoms in the structure and output
def write_pdb_without_H(topology, positions, file_output):
    mod = Modeller(topology, positions)
    hydrogens = [a for a in mod.topology.atoms() if a.element == elem.hydrogen]
    if hydrogens:
        mod.delete(hydrogens)
    with open(file_output, 'w+') as fout:
        PDBFile.writeFile(mod.topology, mod.positions, fout, keepIds=True)

# modify HID/HIP/HIE to HIS
def rename_HIS_variants(file_pdb, enabled=True):
    if not enabled:
        return
    lines = []
    with open(file_pdb, 'r') as fpdb:
        for line in fpdb:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                resname = line[17:20]
                if resname in ("HID", "HIE", "HIP"):
                    line = line[:17] + "HIS" + line[20:]
            lines.append(line)
    with open(file_pdb, 'w') as fout:
        fout.writelines(lines)


############    MAIN    ############
def main():
    if len(sys.argv) != 4:
        print("Usage: python refine_model_openmm.py <input.pdb> <steps> <output.pdb>")
        sys.exit(1)

    input_pdb = sys.argv[1]
    num_steps = max(0, int(sys.argv[2]))
    output_pdb = sys.argv[3]

    pdb_fixer = PDBFixer(filename=input_pdb)
    pdb_fixer.findMissingResidues()
    pdb_fixer.findMissingAtoms()
    pdb_fixer.addMissingAtoms()

    modeller = Modeller(pdb_fixer.topology, pdb_fixer.positions)
    modeller.addHydrogens()

    ff = try_forcefield_with_phos()
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )
    add_backbone_position_restraints(system, modeller.topology, modeller.positions, 1.0)

    md_platform = Platform.getPlatformByName("CPU")
    integrator = VerletIntegrator(0.002 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator, md_platform)
    sim.context.setPositions(modeller.positions)

    tolerance = 10.0 * unit.kilojoule_per_mole / unit.nanometer
    LocalEnergyMinimizer.minimize(sim.context, tolerance, num_steps)

    final_stage = sim.context.getState(getPositions=True)
    final_positions = final_stage.getPositions()
    write_pdb_without_H(modeller.topology, final_positions, output_pdb)
    rename_HIS_variants(output_pdb, True)


if __name__ == "__main__":
    main()
