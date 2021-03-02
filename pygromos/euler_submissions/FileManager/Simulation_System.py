"""
    This module carries the System class
"""
import glob
import os,warnings
from typing import List

from pygromos.euler_submissions.FileManager._base import _fileManagment_base_class
from pygromos.euler_submissions.FileManager import Solvents, Topology
from pygromos.utils.bash import check_path_dependencies, copy_file


class System(_fileManagment_base_class):
    """System

        This class is the representative of a system for a simulation. it should carry all needed requirements for MD simulation.
    """
    name:str = None
    top:Topology = None
    coordinates:List[str] = None
    solvent:Solvents.Solvent = None
    steps = None    #TODO: remvoe? not USED!

    def __init__(self, top: Topology =None, coordinates: List[str]=None, solvent: Solvents.Solvent=None, name:str = ""):
        if(name!=None and top != None and coordinates != None):
            self.name = name
            self.top = top
            self.coordinates = coordinates
            self.solvent = solvent
        else:
            raise IOError("Please give either name or dictionary attribute to "+self.__class__.__name__+".__init__(...)")

    def check(self):
        """Checks if all paths of this system are valid. :return:"""
        check_paths = self.coordinates
        check_paths += self.top._return_all_paths()
        check_paths += self.solvent._return_all_paths()
        check_path_dependencies(check_paths)

    def move_input_coordinates(self, coordinate_folder:str, do_not_copy:bool=False,coord_format:str=".cnf"):
        new_seed_list = []
        if(isinstance(self.coordinates, str) and os.path.isdir(self.coordinates)):
            cnfs=glob.glob(self.coordinates+"/*"+coord_format)
            if(len(cnfs)==0):
               raise IOError("Could not find any .cnf here", self.coordinates)
        elif(isinstance(self.coordinates, str) and os.path.isfile(self.coordinates)):
            cnfs = [self.coordinates]

        else:
            cnfs = self.coordinates


        for cnf in cnfs:
            out_cnf = coordinate_folder + "/" + os.path.basename(cnf)
            if(do_not_copy):
                new_seed_list.append(out_cnf)
            elif(os.path.exists(out_cnf)):
                warnings.warn("did not copy Cnf-file, Path already exists: " + str(out_cnf))
                new_seed_list.append(out_cnf)
            elif(os.path.exists(cnf)):
                new_seed_list.append(copy_file(cnf, out_cnf))
            else:
                raise IOError("S-OPT: Expected that coordinate files exist (do_not_copy:"+str(do_not_copy)+"! \n" + str(
                self.coordinates))

        if(isinstance(self.coordinates, str) and  os.path.isfile(self.coordinates)):
            self.coordinates=out_cnf
        else:
            self.coordinates = new_seed_list
