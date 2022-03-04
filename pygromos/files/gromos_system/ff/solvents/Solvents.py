"""
    This module collects all solvent subclasses for python file and information management with gromos

"""


class Solvent:
    """Solvent
    This class is giving the needed solvent infofmation for gromos in an obj.
    #TODO: CONFs & TOPO in data
    """

    name: str = None
    coord_file_path: str = None

    def __init__(self, name: str = None, coord_file: str = None):

        if name != None:
            self.name = name
            self.coord_file_path = coord_file
        else:
            raise IOError("DID not get correct Constructor arguments in " + self.__class__.name)

    def _return_all_paths(self) -> list:
        coll = []
        if self.coord_file_path != None:
            coll.append(self.coord_file_path)
        return coll


class H2O(Solvent):
    def __init__(self, coord_file_path: str = None):

        if coord_file_path != None:
            super().__init__(name="H2O")
            self.coord_file_path = coord_file_path
            self.atomNum = 3
        else:
            raise IOError("DID not get correct Constructor arguments in " + self.__class__.name)


class CHCL3(Solvent):
    def __init__(self, coord_file_path: str = None):

        if coord_file_path != None:
            super().__init__(name="CHCL3")
            self.coord_file_path = coord_file_path
            self.atomNum = 5
        else:
            raise IOError("DID not get correct Constructor arguments in " + self.__class__.name)


class DMSO(Solvent):
    def __init__(self, coord_file_path: str = None):

        if coord_file_path != None:
            super().__init__(name="DMSO")
            self.coord_file_path = coord_file_path
            self.atomNum = 4
        else:
            raise IOError("DID not get correct Constructor arguments in " + self.__class__.name)
