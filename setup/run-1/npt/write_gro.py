import MDAnalysis as mda

u = mda.Universe("min.pdb", "output.xtc")
u.trajectory[-1]
u.select_atoms("not resname DUM").write("equil.gro")
