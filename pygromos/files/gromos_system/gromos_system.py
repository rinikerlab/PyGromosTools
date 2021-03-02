"""
File:            gromos++ system folder managment
Warnings: this class is WIP!
TODO:REWORK
Description:
    This is a super class for the collection of all gromos files used inside a project
Author: Marc Lehner
"""

#imports
import warnings, glob, importlib

from rdkit import Chem
from rdkit.Chem import AllChem

from pygromos.files.topology import top
from pygromos.files.blocks.topology_blocks import FORCEFIELD
from pygromos.gromos import gromosPP
from pygromos.data.ff import Gromos54A7
from pygromos.data.ff import Gromos2016H66
from pygromos import files
from pygromos.files.topology.top import Top
from pygromos.files.topology.ifp import ifp
from pygromos.files.coord.cnf import Cnf
from pygromos.files.imd import Imd
from pygromos.data import ff
from pygromos.data.ff import *

if(importlib.util.find_spec("openforcefield") != None):
    from openforcefield.topology import Molecule, Topology
    from openforcefield.typing.engines import smirnoff
    from pygromos.files.gromos_system import openforcefield2gromos
    from pygromos.files.gromos_system import serenityff
    has_openff = True
else:
    has_openff = False


class gromos_system():
    def __init__(self, folder:str, name:str, smiles:str = None, rdkitMol:Chem.rdchem.Mol = None, readIn = False, Forcefield = "2016H66", auto_convert=True):
        self.hasData = False
        self.folder = folder
        self.name = name
        self.smiles = smiles
        self.Forcefield = None
        self.mol = Chem.Mol()
        self.top = Top(input=None)
        self.top._orig_file_path = self.folder + "/" + self.name + ".top"
        self.cnf = Cnf(input=None)
        self.cnf._orig_file_path = self.folder + "/" + self.name + ".cnf"
        self.imd = Imd(input=None)
        self.imd._orig_file_path = self.folder + "/" + self.name + ".imd"

        if smiles==None and rdkitMol==None and readIn==False:
            warnings.warn("No data provided to gromos_system\nmanual work needed")
            

        #import molecule from smiles using rdkit
        if smiles:
            self.mol = Chem.MolFromSmiles(smiles)
            Chem.AddHs(self.mol)
            AllChem.EmbedMolecule(self.mol)
            self.hasData=True        

        #import  molecule from RDKit
        if rdkitMol:
            self.mol = rdkitMol
            self.smiles = Chem.MolToSmiles(self.mol)
            self.hasData=True

        #try to read existing files in
        if readIn:
            self.read_files()
            if True: #TODO: make sanity check on import
                self.hasData=True

        #set forcefield
        self.import_Forcefield(Forcefield=Forcefield)

        # automate all conversions if possible
        if auto_convert and self.hasData:
            if Forcefield == "2016H66" or Forcefield == "54A7":
                pass #TODO: make_top()
            if Forcefield == "smirnoff":
                if (has_openff):
                    raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                      "Please install the package for this feature!")
                else:
                    self.top = openforcefield2gromos(Molecule.from_rdkit(self.mol), self.top, forcefield_name=self.Forcefield).convert_return()
            if Forcefield == "serenityff":
                if (has_openff):
                    raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                      "Please install the package for this feature!")
                else:
                    self.serenityff = serenityff()

    
        
    def import_Forcefield(self, Forcefield:str = "2016H66"):
        if Forcefield == "2016H66":
            self.Forcefield = Gromos2016H66.ifp
            self.ifp = ifp(self.Forcefield)
        elif Forcefield == "54A7":
            self.Forcefield = Gromos54A7.ifp
            self.ifp = ifp(self.Forcefield)
        elif Forcefield == "smirnoff":
            if(has_openff):
                raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                  "Please install the package for this feature!")
            else:
                filelist = glob.glob(ff.data_ff_SMIRNOFF + '/*.xml')
                filelist.sort()
                for f in filelist:
                    try:
                        smirnoff.ForceField(f)
                        self.Forcefield = f
                        break
                    except:
                        pass
        elif Forcefield == "serenityff":
            if (has_openff):
                raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                  "Please install the package for this feature!")
            else:
                raise "WIP"
        

    def read_files(self):
        try:
            self.top.read_file()
        except:
            print("no topology found")
        try:
            self.cnf.read_file()
        except:
            print("no conformation found")
        try:
            self.imd.read_file()
        except:
            print("no GROMOS input file found")

    def write_files(self, Top=True, Cnf=True, Imd=True, Mol=False):
        if Top:
            self.top.write(self.top._orig_file_path)
        if Cnf:
            self.cnf.write(self.cnf._orig_file_path)
        if Imd:
            self.imd.write(self.imd._orig_file_path)
        if Mol:
            print(Chem.MolToMolBlock(self.mol),file=open(self.folder + "/" + self.name + ".mol",'w+'))
            
            

    def rdkitImport(self, inputMol: Chem.rdchem.Mol):
        self.mol = inputMol

    def rdkit2Gromos(self):
        # generate top
        self.top.add_block(blocktitle="TITLE", content=["topology generated with PyGromos from RDKit molecule: " + str(Chem.MolToSmiles(self.mol)) + "    " + str(self.name)],verbose=True)
        self.top.add_block(blocktitle="PHYSICALCONSTANTS", content=[""],verbose=True)
        self.top.add_block(blocktitle="TOPVERSION", content=[["2.0"]])
        


