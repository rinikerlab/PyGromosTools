from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system
from pygromos.simulations.approaches.solvation_free_energy_calculation.solvation_free_energy import Solvation_free_energy_calculation


import argparse
import pandas as pd
from rdkit import Chem

parser = argparse.ArgumentParser(description="Setup Solvation Free Energy Calculation")

parser.add_argument('--id',help='id of molecule to calculate',type=int,default=173)
args = parser.parse_args()

data = pd.read_table('/home/kpaul/pygromos_new_fork/PyGromosTools/combined_list_selected.dat')

# Get desired information
smiles = data[data['#No.'] == args.id].SMILES.values[0]
density = data[data['#No.'] == args.id]['Density_[kg/m^3]'].values[0]

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
number_of_atoms = mol.GetNumAtoms()

# Perform simulation
print(smiles,density,number_of_atoms)
#try:
sf = Solvation_free_energy_calculation(input_system=smiles,
                                       work_folder="/home/kpaul/pygromos_new_fork/PyGromosTools/testbench3",
                                       system_name=smiles, forcefield=forcefield_system("openforcefield"),
                                       density=density, num_molecules=512, num_atoms=number_of_atoms, nmpi=1, nomp=1,
                                       subsystem="local", amberscaling=True)

sf.create_liq()
emin_sys, jobID = sf.minimize_liq(gromos_system=sf.groSys_liq,prev_JobID=-1)
sf.nmpi = 8
eq_sys, jobID = sf.eq_liq(gromos_system=emin_sys,prev_JobID=jobID)
ti_sys, jobID = sf.ti_liq(gromos_system=eq_sys,prev_JobID=jobID,n_points=5)
sf.calculate_solvation_free_energy(n_points=5)

'''
   Lambda        avg       rmsd        ee
0    0.00  31.257825   6.455767  0.482510
1    0.25 -14.987395  14.602506  0.543161
2    0.50  31.861698  24.305928  0.411241
3    0.75  55.893268  45.723446  3.582326
4    1.00  30.283032  18.442700  1.352284
Solvation Free Energy: -24.07397894588648 pm -1.5966018745584094

'''

# except:
#     print(args.id,' failed')
