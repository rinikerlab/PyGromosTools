# compute some dihedrals
from pygromos.files.trajectory import trc

traj = trc.Trc(traj_path="menthol.trc")
print(traj.dihedrals(atom_pairs=[[29, 27, 10, 12]]))
print(traj.dihedrals(atom_pairs=[[4, 10, 12, 15]]))
