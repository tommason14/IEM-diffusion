from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import load_file
from parmed.topologyobjects import Atom
from mdtraj.reporters import XTCReporter

# OpenMM parsers can't read virtual sites
gro = load_file("../npt/equil.gro")
top = load_file("../../pack/topol.top")
top.box = gro.box[:]
system = top.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=12 * angstroms,
    constraints=HBonds,
    rigidWater=True,
)

integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
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
simulation.step(100_000_000)

# Write checkpoint to restart from - for NVT run
# To restart:
# with open('checkput.chk', 'rb') as f:
#     simulation.context.loadCheckpoint(f.read())
