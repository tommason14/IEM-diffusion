from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import load_file
from parmed.topologyobjects import Atom
from mdtraj.reporters import XTCReporter

# OpenMM parsers can't read virtual sites
gro = load_file("../../pack/pack.gro")
top = load_file("../../pack/topol.top")
top.box = gro.box[:]
system = top.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=12 * angstroms,
    constraints=HBonds,
    rigidWater=True,
)

barostat = MonteCarloBarostat(1 * bar, 300 * kelvin)
system.addForce(barostat)
# parmed atom object have no residue name/id, so have to search for residues and then find the atoms
# within those residues
polymers = [r for r in top.residues if r.name == "pol"]
polymer_atoms = [a.idx for pol in polymers for a in pol.atoms]
N_RES = len(top.residues)

# Add position restraints to the polymer chains
restraint = HarmonicBondForce()
restraint.setUsesPeriodicBoundaryConditions(True)
system.addForce(restraint)
nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
positions = gro.positions
for i in polymer_atoms:
    j = system.addParticle(0)
    charge = 0.0 * coulomb
    sigma = 0.0 * nanometer
    eps = 0.0 * kilojoules_per_mole
    nonbonded.addParticle(charge, sigma, eps)
    nonbonded.addException(i, j, charge, sigma, eps)
    # add the dummy atoms to a new residue with resname "DUM"
    newatom = Atom(atomic_number=0, name="DU", type="DU")
    top.add_atom(newatom, "DUM", N_RES + 1)
    restraint.addBond(i, j, 0 * nanometers, 100 * kilojoules_per_mole / nanometer ** 2)
    positions.append(positions[i])

integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()
state = simulation.context.getState(getPositions=True)
with open("min.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
simulation.reporters.append(XTCReporter("output.xtc", 5000))
simulation.reporters.append(CheckpointReporter("output.chk", 5000))  # restart
simulation.reporters.append(
    StateDataReporter(
        "output.csv",
        1000,
        step=True,
        time=True,
        kineticEnergy=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        volume=True,
        speed=True,
    )
)
simulation.step(5_000_000)

# Write checkpoint to restart from - for NVT run
# To restart:
# with open('checkput.chk', 'rb') as f:
#     simulation.context.loadCheckpoint(f.read())
