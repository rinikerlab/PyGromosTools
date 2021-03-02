from typing import List, Dict
from pygromos.files._basics import parser, _general_gromos_file


class ene_ana_lib(_general_gromos_file._general_gromos_file):
    """in_ene_ana_lib
        This class represents the ene_ana_lib_file from gromos. (So if you want to understand the structure and use, please check gromos doc)
        It inherits from general_gromos_file.

    """
    required_blocks:List[str] = ["TITLE"]
    block_names:Dict[str, str] = {"TITLE": "title_block", "ENEVERSION":"ene_ana_version",
                                  "ENERTRJ":"ene_traj_structure", "FRENERTRJ": "free_traj_structure", "VARIABLES":"ene_ana_variables" }


    ENETRJ_structure_before_2019 =None

    def __init__(self, in_path: str):
        if(type(in_path) is str):
            self.read_ene_ana_lib_file(in_path)
            self.path = in_path
        elif(type(in_path) is dict):
            self.data = in_path
            self.data = None
        else:
            raise IOError("in_data is none of the expected types str or dict")

    def read_ene_ana_lib_file(self, path):
        print("Begin read in of in_ene_ana_lib-file: "+path)
        data = parser.read_ene_ana_lib(path)
        for key, sub_content in data.items():
            #print(sub_content)
            try:
                self.add_block(blocktitle=key, content=sub_content)
            except Exception as err:
                    raise IOError("could not read in in_ene_ana_lib block!\n"+"\n".join(err.args))
        print("End of read in of in_ene_ana_lib-file")
        return 0

    def check_sanity(self):
        raise NotImplementedError("this is not implemented")

    def convert_gromos_old_2019(self):
        raise NotImplementedError("this is not implemented")


ene_ana_free_Energy_struct_old = (
    "FRENERTRJ"
    "# block definition for the free energy trajectory file."
    "# which is specified by the input flag fr_files of program ene_ana."
    "#"
    "# syntax as for the ENERTRJ definition"
    "#"
    "# Following is the definition for a gromosXX free energy trajectory."
    "#"
    "  block TIMESTEP"
    "    subblock TIME 2 1"
    "  block FREEENERDERIVS03"
    "    subblock RLAM 1 1"
    "    subblock FREEENER 38 1"
    "    size NUM_BATHS"
    "    subblock FREEKINENER NUM_BATHS 1"
    "    size NUM_ENERGY_GROUPS"
    "    subblock FREEBONDED NUM_ENERGY_GROUPS 5"
    "    subblock FREENONBONDED matrix_NUM_ENERGY_GROUPS 4"
    "    subblock FREESPECIAL NUM_ENERGY_GROUPS 9"
    "    size NUM_EDS_STATES"
    "    subblock FREEEDS NUM_EDS_STATES 2 "
    "END"
)
ene_ana_ene_traj_struct_old=(
    "ENERTRJ\n"
    "# block definition for the energy trajectory file.\n"
    "# which is specified by the input flag en_files of program ene_ana.\n"
    "#\n"
    "# Use keyword 'block' to specify the _blocks\n"
    "#             'subblock' to specify name and dimensions of a set of data\n"
    "#             'size' to specify a size that should be read in from the file\n"
    "#                    this size can be used as dimension specification\n"
    "#                    in a subblock definition. Using the prefix 'matrix_'\n"
    "#                    with such a definition will expand the size N to\n"
    "#                    N*(N+1)/2\n"
    "#\n"
    "# Following is the definition for a gromosXX energy trajectory\n"
    "#\n"
    "  block TIMESTEP\n"
    "    subblock TIME 2 1\n"
    "  block ENERGY03\n"
    "    subblock ENER 38 1\n"
    "    size NUM_BATHS\n"
    "    subblock KINENER NUM_BATHS 3\n"
    "    size NUM_ENERGY_GROUPS\n"
    "    subblock BONDED NUM_ENERGY_GROUPS 5\n"
    "    subblock NONBONDED matrix_NUM_ENERGY_GROUPS 4\n"
    "    subblock SPECIAL NUM_ENERGY_GROUPS 11\n"
    "    size NUM_EDS_STATES\n"
    "    subblock EDS NUM_EDS_STATES 3\n"
    "  block VOLUMEPRESSURE03\n"
    "    subblock MASS 1 1\n"
    "    size NUM_BATHS\n"
    "    subblock TEMPERATURE NUM_BATHS 4\n"
    "    subblock VOLUME 10 1\n"
    "END\n"
)