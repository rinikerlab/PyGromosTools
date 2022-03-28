import os
import shutil
import subprocess
from rdkit import Chem


from pygromos.files.coord.cnf import Cnf
from pygromos.files.forcefield._generic_force_field import _generic_force_field
from pygromos.files.topology.top import Top

from pygromos.gromos import GromosPP
from pygromos.data import topology_templates
from pygromos import data

from pygromos.utils.typing import List


class AmberFF(_generic_force_field):
    # static variables to control solvation
    solvate = False
    solventbox = "TIP3PBOX"

    # don't remove temporary directories after execution by default
    clean = False

    def __init__(self, name: str = "amber", path_to_files: str = None, auto_import: bool = True, verbose: bool = False):
        super().__init__(name, path_to_files=path_to_files, auto_import=auto_import, verbose=verbose)
        if auto_import:
            self.auto_import_ff()

    def auto_import_ff(self):
        # check path
        if self.path_to_files is not None:
            if (
                isinstance(self.path_to_files, List)
                and len(self.path_to_files) > 0
                and isinstance(self.path_to_files[0], str)
            ):
                self.amber_basedir = self.path_to_files[0]
        elif shutil.which("tleap") is not None:
            self.amber_basedir = os.path.abspath(os.path.dirname(shutil.which("tleap")) + "/../")
        else:
            raise ImportError(
                "Could not import GAFF FF as ambertools was missing! " "Please install the package for this feature!"
            )

        if self.verbose:
            print("Found amber: " + str(self.amber_basedir))

        self.amber_bindir = self.amber_basedir + "/bin"
        self.leaprc_files = [
            self.amber_basedir + "/dat/leap/cmd/leaprc.gaff",
            self.amber_basedir + "/dat/leap/cmd/leaprc.water.tip3p",
        ]
        self.frcmod_files = [self.amber_basedir + "/dat/leap/parm/frcmod.chcl3"]

        for leaprc in self.leaprc_files:
            if not os.path.isfile(leaprc):
                raise ImportError("could not find ff file " + leaprc)

        for frcmod in self.frcmod_files:
            if not os.path.isfile(frcmod):
                raise ImportError("could not find ff file " + frcmod)

    def create_top(
        self, mol: str, in_top: Top, in_mol2_file: str = None, work_folder: str = None, gromosPP: GromosPP = None
    ) -> Top:
        self.mol = mol
        self.in_mol2_file = in_mol2_file
        self.work_folder = work_folder
        self.gromosPP = gromosPP

        if in_mol2_file is None:
            self.create_mol2()

        self.amber = amber2gromos(
            in_mol2_file=self.in_mol2_file,
            mol=self.mol,
            forcefield=self,
            gromosPP=self.gromosPP,
            work_folder=work_folder,
            solvate=self.solvate,
            solventbox=self.solventbox,
            clean=self.clean,
        )

        if in_top.path is None:
            self.top = Top(self.amber.get_gromos_topology())
        elif isinstance(in_top, Top):
            self.top = in_top + Top(self.amber.get_gromos_topology())
        elif isinstance(in_top, str):
            self.top = Top(in_top) + Top(self.amber.get_gromos_topology())
        else:
            raise TypeError("in_top is of wrong type")

        return self.top

    def create_cnf(
        self, mol: str, in_cnf: Top = None, work_folder: str = None, in_mol2_file: str = None, gromosPP: GromosPP = None
    ) -> Cnf:
        if self.amber is None:
            self.create_mol2()
            self.amber = amber2gromos(
                in_mol2_file=self.in_mol2_file,
                mol=self.mol,
                forcefield=self,
                gromosPP=self.gromosPP,
                work_folder=work_folder,
                solvate=self.solvate,
                solventbox=self.solventbox,
                clean=self.clean,
            )
        if in_cnf.path is None:
            self.cnf = Cnf(self.amber.get_gromos_coordinate_file())
        elif isinstance(in_cnf, Cnf):
            self.cnf = in_cnf + Cnf(self.amber.get_gromos_coordinate_file())
        elif isinstance(in_cnf, str):
            self.cnf = Cnf(in_cnf) + Cnf(self.amber.get_gromos_coordinate_file())
        else:
            raise TypeError("in_cnf is of wrong type")

        return self.cnf

    def create_mol2(self, mol: str = None):
        pass


class amber2gromos:
    def __init__(
        self,
        in_mol2_file: str,
        mol: Chem.rdchem.Mol,
        gromosPP: GromosPP,
        forcefield: _generic_force_field,
        work_folder: str = ".",
        solvate: bool = False,
        solventbox: str = None,
        clean: bool = False,
    ):
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
        solvate: bool, optional
            should the topology be solvated? (default: False)
        solventbox: str, optional
            what solvent should be used for solvation? e.g. TIP3PBOX or CHCL3BOX (default: TIP3PBOX)
        clean: bool, optional
            should temporary ambertool files be removed after parameterization? (default: False)
        """

        self.in_mol2_file = os.path.abspath(in_mol2_file)

        current_dir = os.getcwd()
        os.chdir(work_folder)

        self.mol = mol
        self.Forcefield = forcefield
        self.gromosPP = gromosPP
        self.solvate = solvate
        self.solventbox = solventbox

        self.mol2_name = "".join(os.path.basename(self.in_mol2_file).split(".")[:-1])

        # ambertools pipeline to generate pdb, crd and prm files
        self.antechamber()
        self.parmchk()
        self.tleap()

        # convert AMBER files to GROMOS
        self.amber2gromos()

        if clean:
            self.cleanup

        os.chdir(current_dir)

    def antechamber(self):
        """
        executes the ambertools program antechamber
        """
        self.antechamber_dir = "antechamber_tmp"
        os.makedirs(self.antechamber_dir, exist_ok=True)
        self.antechamber_dir = os.path.abspath(self.antechamber_dir)
        os.chdir(self.antechamber_dir)

        net_charge = Chem.rdmolops.GetFormalCharge(self.mol)

        self.antechamber_out_mol2 = self.mol2_name + ".mol2"

        command = [self.Forcefield.amber_bindir + "/antechamber"]
        command.extend(["-i", self.in_mol2_file])
        command.extend(["-fi", "mol2"])
        command.extend(["-o", self.antechamber_out_mol2])
        command.extend(["-fo", "mol2"])
        command.extend(["-s", "2"])
        command.extend(["-c", "bcc"])
        command.extend(["-nc", str(net_charge)])

        print(command)
        success = not subprocess.call(command)
        if not success:
            raise RuntimeError("execution of parmchk2 not successful for: " + " ".join(command))

        os.chdir("..")

    def parmchk(self):
        """
        executes the ambertools program parmchk
        """
        self.parmchk_dir = "parmchk_tmp"
        os.makedirs(self.parmchk_dir, exist_ok=True)
        self.parmchk_dir = os.path.abspath(self.parmchk_dir)
        os.chdir(self.parmchk_dir)

        self.parmchk_out_frcmod = self.mol2_name + ".frcmod"

        command = [self.Forcefield.amber_bindir + "/parmchk2"]
        command.extend(["-i", self.antechamber_dir + "/" + self.antechamber_out_mol2])
        command.extend(["-f", "mol2"])
        command.extend(["-o", self.parmchk_out_frcmod])

        print(command)
        success = not subprocess.call(command)
        if not success:
            raise RuntimeError("execution of parmchk2 not successful for: " + " ".join(command))
        os.chdir("..")

    def tleap(self):
        """
        executes the ambertools program tleap
        """
        self.tleap_dir = "tleap_tmp"
        os.makedirs(self.tleap_dir, exist_ok=True)
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
        cmd_file.write(
            "frcmod_" + self.mol2_name + " = loadamberparams " + self.parmchk_dir + "/" + self.parmchk_out_frcmod + "\n"
        )
        cmd_file.write(self.mol2_name + " = loadmol2 " + self.antechamber_dir + "/" + self.antechamber_out_mol2 + "\n")

        if self.solvate:
            cmd_file.write("solvateBox " + self.mol2_name + " " + self.solventbox + " 14 0.75\n")
            self.prm_file = "".join(self.prm_file.split(".")[:-1]) + "_" + self.solventbox + ".prm"
            self.crd_file = "".join(self.crd_file.split(".")[:-1]) + "_" + self.solventbox + ".crd"
            self.pdb_file = "".join(self.pdb_file.split(".")[:-1]) + "_" + self.solventbox + ".pdb"

        cmd_file.write("saveamberparm " + self.mol2_name + " " + self.prm_file + " " + self.crd_file + "\n")
        cmd_file.write("savepdb " + self.mol2_name + " " + self.pdb_file + "\n")
        cmd_file.write("quit\n")
        cmd_file.close()

        command = [self.Forcefield.amber_bindir + "/tleap"]
        command.extend(["-f", cmd_file_name])

        print(command)
        success = not subprocess.call(command)
        if not success:
            raise RuntimeError("execution of parmchk2 not successful for: " + " ".join(command))
        os.chdir("..")

    def amber2gromos(self):
        """
        converts an AMBER (GAFF) parameter and crd file to a GROMOS top and cnf file
        """
        if self.solvate:
            self.gromos_topology = self.mol2_name + "_" + self.solventbox + ".top"
        else:
            self.gromos_topology = self.mol2_name + ".top"
        spc_template = os.path.dirname(os.path.abspath(topology_templates.__file__)) + "/spc.top"
        self.gromosPP.amber2gromos(ambertop=self.prm_file, solvent=spc_template, out_path=self.gromos_topology)
        self.gromos_topology = os.path.abspath(self.gromos_topology)
        print("converted topology saved to " + self.gromos_topology)

        pdb2g96_lib = os.path.dirname(os.path.abspath(data.__file__)) + "/pdb2g96.lib"
        if self.solvate:
            self.gromos_coordinate_file = self.mol2_name + "_" + self.solventbox + ".cnf"
        else:
            self.gromos_coordinate_file = self.mol2_name + ".cnf"
        self.gromosPP.pdb2gromos(
            in_top_path=self.gromos_topology,
            in_pdb_path=self.pdb_file,
            in_lib_path=pdb2g96_lib,
            out_cnf_path=self.gromos_coordinate_file,
        )

        self.gromos_coordinate_file = os.path.abspath(self.gromos_coordinate_file)

        if self.solvate:
            with open(self.crd_file, "r") as f:
                last_line = f.readlines()[-1]
            print(last_line)
            x = float(last_line.split()[0]) * 0.1  # x-coordinate in nm
            y = float(last_line.split()[1]) * 0.1  # y-coordinate in nm
            z = float(last_line.split()[2]) * 0.1  # z-coordinate in nm

            a1 = last_line.split()[3]  # box angle
            a2 = last_line.split()[4]  # box angle
            a3 = last_line.split()[5]  # box angle
            shutil.copyfile(self.gromos_coordinate_file, "tmp_cnf")
            file_out = open(self.gromos_coordinate_file, "w")
            file_in = open("tmp_cnf")

            for line in file_in:
                if "GENBOX" not in line:
                    file_out.write(line)
                else:
                    file_out.write("GENBOX\n")
                    file_out.write("    1\n")
                    file_out.write(
                        "    " + str(float(x) / 10) + "  " + str(float(y) / 10) + "  " + str(float(z) / 10) + "\n"
                    )
                    file_out.write("    " + str(a1) + "  " + str(a2) + "  " + str(a3) + "\n")
                    file_out.write("    0.000000000    0.000000000    0.000000000\n")
                    file_out.write("    0.000000000    0.000000000    0.000000000\n")
                    file_out.write("END")
                    break

            file_in.close()
            file_out.close()
            os.remove("tmp_cnf")

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
        return self.gromos_coordinate_file

    def cleanup(self):
        """
        removes temporary parmchk, antechamber, and tleap directories
        """
        os.rmdir(self.tleap_dir)
        os.rmdir(self.antechamber_dir)
        os.rmdir(self.parmchk_dir)
