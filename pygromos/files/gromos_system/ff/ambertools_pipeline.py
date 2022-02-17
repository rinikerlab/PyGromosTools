"""
File:            ambertools_pipeline
Warnings: this class is WIP!
Description:
    This is a class to parameterize a molecule with ambertools
Author: Salomé Rieder
"""

#imports
import importlib
from pygromos.files.topology.top import Top
from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system

from pygromos.data import topology_templates
import os
from rdkit import Chem

class ambertools_pipeline():
    def __init__(self, in_mol2_file: str, mol: Chem.rdchem.Mol, forcefield: forcefield_system = None):

        self.in_mol2_file = in_mol2_file
        self.mol = mol
        self.Forcefield = forcefield

        self.antechamber()
        self.parmchk()
        self.tleap()

    def antechamber(self):
        self.antechamber_dir = "antechamber_tmp"
        os.makedirs(self.antechamber_dir, exist_ok = True)
        self.antechamber_dir = os.path.abspath(self.antechamber_dir)
        os.chdir(self.antechamber_dir)

        net_charge = Chem.rdmolops.GetFormalCharge(self.mol)
        self.antechamber_mol2 = os.path.basename(self.in_mol2_file)

        command = self.Forcefield.amber_bindir + "/antechamber"\
                    + " -i " + self.in_mol2_file \
                    + " -fi mol2" \
                    + " -o " + self.antechamber_mol2 \
                    + " -fo mol2 -s 2 -c bcc " \
                    + " -nc " + str(net_charge)

        print(command)
        os.system(command)
        os.chdir('..')

    def parmchk(self):
        self.parmchk_dir = "parmchk_tmp"
        os.makedirs(self.parmchk_dir, exist_ok = True)
        self.parmchk_dir = os.path.abspath(self.parmchk_dir)
        os.chdir(self.parmchk_dir)

        self.parmchk_frcmod = ''.join(self.antechamber_mol2.split('.')[:-1]) + ".frcmod"

        command = self.Forcefield.amber_bindir + "/parmchk2" \
                    + " -i " + self.antechamber_dir \
                    + "/" + os.path.basename(self.in_mol2_file) \
                    + " -f mol2 -o " + self.parmchk_frcmod

        print(command)
        os.system(command)
        os.chdir('..')

    def tleap(self):
        self.tleap_dir = "tleap_tmp"
        os.makedirs(self.tleap_dir, exist_ok = True)
        self.tleap_dir = os.path.abspath(self.tleap_dir)
        os.chdir(self.tleap_dir)

        cmd_file = open("tleap.cmd", "w")
        for leaprc in self.Forcefield.leaprc_files:
            cmd_file.write("source " + leaprc + "\n")

        cmd_file.write("set default pbradii mbondi\n")
        cmd_file.write("set default nocenter on\n")

        #command = self.Forcefield.amber_bindir + "/tleap"

        #print(command)
        #os.system(command)
        os.chdir('..')
