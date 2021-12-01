"""
File:            forcefield managment for gromos_system
Warnings: this class is WIP!

Description:
    This class handles all possible forcefields gromos_system accepts
Author: Marc Lehner
"""

# imports
import glob
import importlib
import collections

from pygromos.files.topology.ifp import ifp
from pygromos.files.topology.top import Top

from pygromos.data import topology_templates
from pygromos.data import ff
from pygromos.data.ff import Gromos2016H66
from pygromos.data.ff import Gromos54A7

if (importlib.util.find_spec("openff") != None):
    from openff.toolkit.typing.engines import smirnoff
    has_openff = True
else:
    has_openff = False



class forcefield_system():
    def __init__(self, name:str="2016H66", path:str=None, auto_import:bool = True):
        self.name = name
        self.path = path
        self.mol_name=None
        if auto_import:
            self.auto_import_ff()

    def auto_import_ff(self):
        if self.name == "2016H66":
            self.path = Gromos2016H66.ifp
            self.ifp = ifp(self.path)
            self.mtb_path = Gromos2016H66.mtb
            self.mtb_orga_path = Gromos2016H66.mtb_orga

        elif self.name == "54A7":
            self.path = Gromos54A7.ifp
            self.ifp = ifp(self.path)
            self.mtb_path = Gromos54A7.mtb

        elif self.name == "openforcefield" or self.name == "smirnoff" or self.name == "off":
            self.import_off()

        elif self.name == "serenityff":
            self.import_off()
            self.top = Top(in_value=topology_templates.topology_template_dir + "/blank_template+spc.top")
            self.develop = False
            self.C12_input={}
            self.partial_charges = collections.defaultdict(float)

    def import_off(self):
        if not has_openff:
            raise ImportError("Could not import smirnoff FF as openFF toolkit was missing! "
                                "Please install the package for this feature!")
        if self.path != None:
            try:
                self.off = smirnoff.ForceField(self.path)
            except:
                raise ImportError("Could not import a OpenForceField from path: " +str(self.path))
        else:
            filelist = glob.glob(ff.data_ff_SMIRNOFF + '/*.offxml')
            filelist.sort()
            filelist.reverse()
            for f in filelist:
                try:
                    self.off = smirnoff.ForceField(f)
                    self.path = f
                    break
                except Exception as err:
                    pass
        print("Found off: "+str(self.path))
