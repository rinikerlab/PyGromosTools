from typing import Union, List, Dict
from pygromos.files.blocks import _general_blocks
from collections import defaultdict


##ResidueLibBlocks

class RESIDUENAMELIB(_general_blocks._generic_gromos_block):
    fields = ["pdb_name", "g96top_name"]
    pdb_top:defaultdict(list)

    def __init__(self, data:Union[List[str], Dict[str, List[str]]]):
        super().__init__(self.__class__.__name__, used=True)

        #print(data)
        if(isinstance(data, List)):
            self.pdb_top= defaultdict(list)
            tuples =  list(map(lambda x: x.strip().split(), filter(lambda x: not x.startswith("#"), data)))
            for pdb, top in tuples:
                self.pdb_top[pdb].append(top)

        elif(isinstance(data, Dict)):
            self.pdb_top = data
        else:
            raise ValueError("did not get appropriate data type.")

        #print(self.pdb_top)

    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result+= self.comment_char+self.field_seperator.join(self.fields)+self.line_seperator

        for pdb, top in sorted(self.pdb_top.items()):
            for topI in top:
                result += self.field_seperator+pdb+self.field_seperator+topI+self.line_seperator
        result += "END" + self.line_seperator
        return result

class ATOMNAMELIB(_general_blocks._generic_gromos_block):
    fields = ["top_resn", "atom_pdb", "atom_top"]
    pdb_top:defaultdict(dict)

    def __init__(self,  data:Union[List[str], Dict[str, Dict[str, str]]]):
        super().__init__(self.__class__.__name__, used=True)

        #print(data)
        if(isinstance(data, List)):
            tuples = [(resn, pdb, top) for resn, pdb, top in list(map(lambda x: x.strip().split(), filter(lambda x: not x.startswith("#"), data)))]
            pdb_top = defaultdict(dict)
            for resn, pdb, top in tuples:
                pdb_top[resn].update({pdb: top})
            self.pdb_top = pdb_top
        elif(isinstance(data, Dict)):
            self.pdb_top = data
        else:
            raise ValueError("did not get appropriate data type.")
        #print(self.pdb_top)


    def block_to_string(self) -> str:
        result = self.name + self.line_seperator
        result+= self.comment_char+self.field_seperator.join(self.fields)+self.line_seperator

        for resn in self.pdb_top:
            atom_dict = self.pdb_top[resn]
            for pdb, top in sorted(atom_dict.items()):
                result += self.field_seperator+resn+self.field_seperator+pdb+self.field_seperator+top+self.line_seperator
        result += "END" + self.line_seperator
        return result