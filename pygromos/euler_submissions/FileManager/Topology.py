"""
    This module contains the topology filemanagment class.
        Todo: maybe rename.
"""

from pygromos.euler_submissions.FileManager._base import _fileManagment_base_class
from pygromos.utils.bash import check_path_dependencies


class Topology(_fileManagment_base_class):
    """Topology
        This is a class collecting all topology files for gromos.
            So you only need to pass one obj. instead of multiple file-paths.
    """
    top_path = None
    disres_path = None
    pertubation_path = None
    def __init__(self, top_path:str=None, disres_path:str=None, pertubation_path:str=None):
        if (top_path != None):
            self.top_path = top_path
            self.disres_path = disres_path
            self.pertubation_path = pertubation_path
        else:
            raise IOError("DID not get correct Constructor arguments in "+self.__class__.name)



    def _return_all_paths(self)->list:
        coll = []
        if(self.top_path != None):
            coll.append(self.top_path)
        if(self.pertubation_path != None):
            coll.append(self.pertubation_path)
        if (self.disres_path != None):
            coll.append(self.disres_path)

        return coll

    def check(self):
        check_path_dependencies(self._return_all_paths())