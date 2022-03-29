from pygromos.files._basics import _general_gromos_file, parser
from pygromos.utils.typing import Union, Dict


class Distance_restraints(_general_gromos_file._general_gromos_file):
    required_blocks = ["TITLE", "DISTANCERESPEC"]
    _gromos_file_ending: str = "disres"

    def __init__(self, in_value: Union[str, Dict] = None):
        self.blocksset = []
        self.block_names = {"TITLE": "title_block", "DISTANCERESSPEC": "distance_res_spec_block"}
        super().__init__(in_value=in_value)

        """
        if(type(path) is str):
            self.path = path
            self.read_disres(path)

        elif(path==None):
            warnings.warn("Warning!: generated empty disres obj!")
        else:
            raise IOError("disres class got "+str(type(path))+" as input. Unknown input type for disres.")
        """

    def read_blocks(self):
        # parse file into dicts
        data = parser.read_disres(self.path)
        # add _blocks as attribute to objects
        for key, sub_content in data.items():
            # print(sub_content)
            self.add_block(block=sub_content)


class Disres(Distance_restraints):
    pass
