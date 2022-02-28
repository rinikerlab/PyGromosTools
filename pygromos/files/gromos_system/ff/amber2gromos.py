"""
File:            amber2gromos
Description:
    This is a class to parameterize a molecule with ambertools. Generates a GROMOS top and cnf file.
Author: Salomé Rieder
"""

#imports
import subprocess
from pygromos.files.topology.top import Top
from pygromos.files.gromos_system.ff.forcefield_system import forcefield_system

from pygromos.data import topology_templates
from pygromos import data
from pygromos.gromos import GromosPP
import os
from rdkit import Chem

class amber2gromos():
    # static variables to control solvation
    solvate = False
    solventbox = "TIP3PBOX"

    # don't remove temporary directories after execution by default
    clean = False

    def __init__(self, 
                    in_mol2_file: str, 
                    mol: Chem.rdchem.Mol, 
                    gromosPP: GromosPP,
                    forcefield: forcefield_system,
                    work_folder: str = "."):
        """
        uses the ambertools programs antechamber, parmchk, and tleap together with
        the GROMOS++ program amber2gromos to generate a GROMOS topology and coordinate file
        for a given molecule

        Parameters
        ----------
        in_mol2_file : str
            mol2 file for molecule to be parameterized
        mol: Chem.rdchem.Mol
            rdkit molecule of the molecule in in_mol2_file
        gromosPP: GromosPP
        forcefield: forcefield_system
        work_folder: str, optional
            where to generate the topology + cnf (default: ".")
        """

        self.in_mol2_file = os.path.abspath(in_mol2_file)

        current_dir = os.getcwd()
        os.chdir(work_folder)

        self.mol = mol
        self.Forcefield = forcefield
        self.gromosPP = gromosPP

        self.mol2_name = ''.join(os.path.basename(self.in_mol2_file).split('.')[:-1])

        # ambertools pipeline to generate pdb, crd and prm files
        self.antechamber()
        self.parmchk()
        self.tleap()

        # convert AMBER files to GROMOS
        self.amber2gromos()

        if(self.clean):
            self.cleanup

        os.chdir(current_dir)

    def antechamber(self):
        """
        executes the ambertools program antechamber
        """
        self.antechamber_dir = "antechamber_tmp"
        os.makedirs(self.antechamber_dir, exist_ok = True)
        self.antechamber_dir = os.path.abspath(self.antechamber_dir)
        os.chdir(self.antechamber_dir)

        net_charge = Chem.rdmolops.GetFormalCharge(self.mol)
        
        self.antechamber_out_mol2 = self.mol2_name + ".mol2"

        command = [self.Forcefield.amber_bindir + "/antechamber"]
        command.extend(["-i",self.in_mol2_file])
        command.extend(["-fi","mol2"])
        command.extend(["-o",self.antechamber_out_mol2])
        command.extend(["-fo","mol2"])
        command.extend(["-s","2"])
        command.extend(["-c","bcc"])
        command.extend(["-nc",str(net_charge)])

        print(command)
        success = not subprocess.call(command)
        if(not success):
            raise RuntimeError("execution of parmchk2 not successful for: " + ' '.join(command))

        os.chdir('..')

    def parmchk(self):
        """
        executes the ambertools program parmchk
        """
        self.parmchk_dir = "parmchk_tmp"
        os.makedirs(self.parmchk_dir, exist_ok = True)
        self.parmchk_dir = os.path.abspath(self.parmchk_dir)
        os.chdir(self.parmchk_dir)

        self.parmchk_out_frcmod = self.mol2_name + ".frcmod"

        command = [self.Forcefield.amber_bindir + "/parmchk2"]
        command.extend(["-i",self.antechamber_dir + "/" + self.antechamber_out_mol2])
        command.extend(["-f","mol2"])
        command.extend(["-o",self.parmchk_out_frcmod])

        print(command)
        success = not subprocess.call(command)
        if(not success):
            raise RuntimeError("execution of parmchk2 not successful for: " + ' '.join(command))
        os.chdir('..')

    def tleap(self):
        """
        executes the ambertools program tleap
        """
        self.tleap_dir = "tleap_tmp"
        os.makedirs(self.tleap_dir, exist_ok = True)
        self.tleap_dir = os.path.abspath(self.tleap_dir)
        os.chdir(self.tleap_dir)

        cmd_file_name = self.mol2_name + ".cmd"
        cmd_file = open(cmd_file_name, "w")
        self.prm_file = self.tleap_dir + "/" + self.mol2_name + ".leap.prm"
        self.crd_file = self.tleap_dir + "/" + self.mol2_name + ".leap.crd"
        self.pdb_file = self.tleap_dir + "/" + self.mol2_name + ".leap.pdb"

        for leaprc in self.Forcefield.leaprc_files:
            cmd_file.write("source " + leaprc + "\n")

        for frcmod in self.Forcefield.frcmod_files:
            cmd_file.write("loadamberparams " + frcmod + "\n")

        cmd_file.write("set default pbradii mbondi\n")
        cmd_file.write("set default nocenter on\n")
        cmd_file.write("frcmod_" + self.mol2_name + " = loadamberparams " + self.parmchk_dir + "/" + self.parmchk_out_frcmod + "\n")
        cmd_file.write(self.mol2_name + " = loadmol2 " + self.antechamber_dir + "/" + self.antechamber_out_mol2 + "\n")

        if(self.solvate):
            cmd_file.write("solvateBox " + self.mol2_name + " " + self.solventbox + " 14 0.75\n")
            self.prm_file = ''.join(self.prm_file.split('.')[:-1]) + "_" + self.solventbox + ".prm"
            self.crd_file = ''.join(self.crd_file.split('.')[:-1]) + "_" + self.solventbox + ".crd"
            self.pdb_file = ''.join(self.pdb_file.split('.')[:-1]) + "_" + self.solventbox + ".pdb"
        
        cmd_file.write("saveamberparm " + self.mol2_name + " " + self.prm_file + " " + self.crd_file + "\n")
        cmd_file.write("savepdb " + self.mol2_name + " " + self.pdb_file + "\n")
        cmd_file.write("quit\n")
        cmd_file.close()

        command = [self.Forcefield.amber_bindir + "/tleap"]
        command.extend(['-f',cmd_file_name])

        print(command)
        success = not subprocess.call(command)
        if(not success):
            raise RuntimeError("execution of parmchk2 not successful for: " + ' '.join(command))
        os.chdir('..')

    def amber2gromos(self):
        """
        converts an AMBER (GAFF) parameter and crd file to a GROMOS top and cnf file
        """
        if(self.solvate):
            self.gromos_topology = self.mol2_name + "_" + self.solventbox + ".top"
        else:
            self.gromos_topology = self.mol2_name + ".top"
        spc_template = os.path.dirname(os.path.abspath(topology_templates.__file__)) + "/spc.top"
        self.gromosPP.amber2gromos(ambertop = self.prm_file, solvent = spc_template, out_path = self.gromos_topology)
        self.gromos_topology = os.path.abspath(self.gromos_topology)
        print("converted topology saved to " + self.gromos_topology)

        pdb2g96_lib = os.path.dirname(os.path.abspath(data.__file__)) + "/pdb2g96.lib"
        if(self.solvate):
            self.gromos_coordinate_file = self.mol2_name + "_" + self.solventbox + ".cnf"
        else:
            self.gromos_coordinate_file = self.mol2_name + ".cnf"
        self.gromosPP.pdb2gromos(in_top_path = self.gromos_topology, in_pdb_path = self.pdb_file, in_lib_path = pdb2g96_lib, out_cnf_path = self.gromos_coordinate_file)
        self.gromos_coordinate_file = os.path.abspath(self.gromos_coordinate_file)
        print("converted coordinates saved to " + self.gromos_coordinate_file)

    def get_gromos_topology(self) -> str:
        """
        Returns
        -------
        str
            returns the path to the converted GROMOS topology
        """
        return self.gromos_topology

    def get_gromos_coordinate_file(self):
        """
        Returns
        -------
        str
            returns the path to the converted GROMOS coordinate file
        """
        print(self.gromos_coordinate_file)
        return self.gromos_coordinate_file

    def cleanup(self):
        """
        removes temporary parmchk, antechamber, and tleap directories
        """
        os.rmdir(self.tleap_dir)
        os.rmdir(self.antechamber_dir)
        os.rmdir(self.parmchk_dir)
