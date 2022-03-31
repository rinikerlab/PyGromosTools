from collections import defaultdict

from pygromos.files.blocks import _general_blocks
from pygromos.utils.typing import Union, List, Dict

# ResidueLibBlocks


class RESIDUENAMELIB(_general_blocks._generic_gromos_block):
    fields = ["pdb_name", "g96top_name"]
    pdb_top: defaultdict(list)

    def __init__(self, content: Union[str, List[str], Dict[str, Dict[str, str]]]):
        """
        content is the content for the  residuenamelib

        Parameters
        ----------
        content : Union[List[str], Dict[str, Dict[str, str]]]
        """
        print(type(content))
        if isinstance(content, list):
            super().__init__(self.__class__.__name__, used=True, content=content)
        elif isinstance(content, str):
            print(content.split())
            super().__init__(self.__class__.__name__, used=True, content=content.split("\n"))
        elif isinstance(content, (dict, defaultdict)):
            super().__init__(self.__class__.__name__, used=True)
            self.pdb_top = content
        else:
            raise ValueError(
                "Class<"
                + str(self.__class__.__name__)
                + ">:\t did not get appropriate data type.\n got: "
                + str(type(content))
                + "\n"
                + str(content)
            )

    def read_content_from_str(self, content: List[str]):
        self.pdb_top = defaultdict(list)
        tuples = list(map(lambda x: x.strip().split(), filter(lambda x: not x.startswith("#"), content)))
        for pdb, top in tuples:
            self.pdb_top[pdb].append(top)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += self.comment_char + self.field_seperator.join(self.fields) + self.line_seperator

        for pdb, top in sorted(self.pdb_top.items()):
            for topI in top:
                result += self.field_seperator + pdb + self.field_seperator + topI + self.line_seperator
        result += "END" + self.line_seperator
        return result


class ATOMNAMELIB(_general_blocks._generic_gromos_block):
    fields = ["top_resn", "atom_pdb", "atom_top"]
    pdb_top: defaultdict(dict)

    def __init__(self, content: Union[List[str], Dict[str, Dict[str, str]]]):
        """
            atom name translation lib part of the resn lib

        Parameters
        ----------
        content : Union[List[str], Dict[str, Dict[str, str]]]
        """
        print(type(content))
        if isinstance(content, list):
            super().__init__(self.__class__.__name__, used=True, content=content)
        elif isinstance(content, str):
            print(content.split())
            super().__init__(self.__class__.__name__, used=True, content=content.split("\n"))
        elif isinstance(content, (dict, defaultdict)):
            super().__init__(self.__class__.__name__, used=True)
            self.pdb_top = content
        else:
            raise ValueError(
                "Class<"
                + str(self.__class__.__name__)
                + ">:\tdid not get appropriate data type.\n got: "
                + str(type(content))
                + "\n"
                + str(content)
            )

    def read_content_from_str(self, content: Union[List[str], Dict[str, Dict[str, str]]]):
        self.pdb_top = defaultdict(list)
        tuples = [
            (resn, pdb, top)
            for resn, pdb, top in list(
                map(lambda x: x.strip().split(), filter(lambda x: not x.startswith("#"), content))
            )
        ]
        pdb_top = defaultdict(dict)
        for resn, pdb, top in tuples:
            pdb_top[resn].update({pdb: top})
        self.pdb_top = pdb_top

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result += self.comment_char + self.field_seperator.join(self.fields) + self.line_seperator

        for resn in self.pdb_top:
            atom_dict = self.pdb_top[resn]
            for pdb, top in sorted(atom_dict.items()):
                result += (
                    self.field_seperator
                    + resn
                    + self.field_seperator
                    + pdb
                    + self.field_seperator
                    + top
                    + self.line_seperator
                )
        result += "END" + self.line_seperator
        return result
