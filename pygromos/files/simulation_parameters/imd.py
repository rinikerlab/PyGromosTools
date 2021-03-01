"""
FUNCTIONLIB:            gromos++ input file functions
Description:
    in this lib, gromosXX input file mainpulating functions are gathered
Author: Kay Schaller & Benjamin Schroeder
"""
import numpy as np
from numbers import Number
from typing import List, Dict, NamedTuple, Iterable
from collections import namedtuple

from pygromos.files._basics import _general_gromos_file, parser
from pygromos.files.blocks import imd_blocks as blocks
from pygromos.utils import amino_acids as aa

def nice_s_vals(svals:Iterable, base10=False) ->list:
    """

    Parameters
    ----------
    svals :
    base10 :

    Returns
    -------

    """

    nicer_labels = []
    if(base10):
        for val in svals:
            if(float(np.log10(val)).is_integer() or val == min(svals)):
                nicer_labels.append(round(val, str(val).count("0")+3))
            else:
                nicer_labels.append("")
    else:
        for val in svals:
            nicer_labels.append(round(float(val), str(val).count("0")+2))
    return nicer_labels

class Imd(_general_gromos_file._general_gromos_file):
    gromos_file_ending:str = "imd"

    _orig_file_path:str
    path:str
    _required_blocks = ["TITLE", "SYSTEM", "STEP", "BOUNDCOND", "FORCE", "CONSTRAINT", "PAIRLIST", "NONBONDED"] # Not jet implemented

    #POSSIBLE GROMOS BLOCKS
    TITLE: blocks.TITLE

    STEP: blocks.STEP
    WRITETRAJ: blocks.WRITETRAJ
    PRINTOUT: blocks.PRINTOUT

    SYSTEM: blocks.SYSTEM
    FORCE: blocks.FORCE

    ENERGYMIN: blocks.ENERGYMIN = None
    EDS:blocks.EDS=None

    CONSTRAINT: blocks.CONSTRAINT = None
    BOUNDCOND: blocks.BOUNDCOND = None
    COMTRANSROT: blocks.COMTRANSROT = None

    PRESSURESCALE: blocks.PRESSURESCALE = None
    MULTIBATH: blocks.MULTIBATH = None

    PAIRLIST: blocks.PAIRLIST = None
    NONBONDED: blocks.NONBONDED = None

    POSITIONRES: blocks.POSITIONRES = None
    DISTANCERES: blocks.DISTANCERES = None
    INITIALISE: blocks.INITIALISE = None

    REPLICA_EDS: blocks.REPLICA_EDS = None
    OLD_REPLICA_EDS: blocks.OLD_REPLICA_EDS = None
    NEW_REPLICA_EDS: blocks.NEW_REPLICA_EDS = None

    def __init__(self, in_value:str, _future_file:bool=False):
        super().__init__(in_value=in_value, _future_file=_future_file)

    def __str__(self):
        text = ""
        if(hasattr(self, "TITLE")):
            text += self.__getattribute__("TITLE").block_to_string()
        for block in sorted(self.get_block_names()):
            if(block == "TITLE" or isinstance(block, type(None))):
                continue
            text += str(self.__getattribute__(block))
        return text

    def read_file(self):
        data = parser.read_imd(self._orig_file_path)
        for key, sub_content in data.items():
            try:
                self.add_block(blocktitle=key, content=sub_content)
            except Exception as err:
                try:
                    print("THE NEW REEDS BLOCK?")
                    if (key == "REPLICA_EDS"):  # TODO: remove as soon as new block is established! or change to old >)
                        self.add_block(blocktitle="NEW_REPLICA_EDS", content=sub_content)
                        self.REPLICA_EDS = self.NEW_REPLICA_EDS
                    else:
                        raise IOError("Could not read in imd " + key + " block!\n values: \n\t" + str(sub_content) + "\n\n" + "\n\t".join(err.args))
                except Exception as err:
                    print("THE OLD REEDS BLOCK?")
                    try:
                        if(key == "REPLICA_EDS"): #TODO: remove as soon as new block is established! or change to old >)
                            self.add_block(blocktitle="OLD_REPLICA_EDS", content=sub_content)
                            self.REPLICA_EDS = self.OLD_REPLICA_EDS
                        else:
                            raise IOError("Could not read in imd "+key+" block!\n values: \n\t"+str(sub_content)+"\n\n"+"\n\t".join(err.args))
                    except Exception as err:
                        raise IOError("could not read in reeds_imd "+key+" block!\n values: \n\t"+str(sub_content)+"\n\n"+"\n\t".join(err.args))
        return {}

    def edit_EDS(self, NUMSTATES:int, S:float, EIR:list, EDS:int=1, ALPHLJ:float=0.0, ALPHC:float=0.0, FUNCTIONAL:int=1):
        if(hasattr(self, "eds_block")):
            self.EDS.ALPHLJ = ALPHLJ
            self.EDS.ALPHC = ALPHC
            self.EDS.FUNCTIONAL = FUNCTIONAL
            self.EDS.NUMSTATES = NUMSTATES
            self.EDS.S = S
            self.EDS.EIR = EIR

        else:
            print("Setting new EDS_block")
            if(type(EIR) == float or type(EIR) == str or type(EIR) == int):
                EIR = [float(EIR) for x in range(NUMSTATES)]

            gromos_name = "EDS"
            setattr(self, gromos_name, blocks.EDS(NUMSTATES, S, EIR, EDS, ALPHLJ, ALPHC, FUNCTIONAL))


    def randomize_seed(self):
        self.INITIALISE.IG = np.random.randint(low=0, high=999999)

    def edit_REEDS(self, REEDS:(bool or int)=None, NUMSTATES:int=None, SVALS: (Number, List[Number])=None, EIR:(Number or Iterable[Number])=None,
                   NRETRIAL:int=None, NREQUIL:int=None, CONT:(bool, int)=None, EDS_STAT_OUT:(bool, int)=None,
                   RETS:List[float]= None, RET:int=None, NATOM:int=None) :   #TODO: old params - to be REMOVED!

        # specific relations are rescued here
        reeds_block = self.REPLICA_EDS
        print(type(reeds_block))

        if(isinstance(reeds_block, blocks.OLD_REPLICA_EDS)):    #old deappriciated-#todo delete if not needed anymore!
            print("GOIN G THROUG OLD IF!")
            if (isinstance(NATOM, (Number, str))):
                reeds_block.NATOM = NATOM
            if (isinstance(RET, (Number, str))):
                reeds_block.RET = RET
            if (isinstance(NUMSTATES, (Number, str)) ):
                reeds_block.NUMSTATES = NUMSTATES

            if (isinstance(SVALS, Iterable)): #edit SVALS
                SVALS = nice_s_vals(SVALS)
                print("SVALS: ", len(SVALS))

                if not isinstance(RETS, Iterable):    #expand to new s value ammount
                    std_RETS_val = reeds_block.RETS[0]
                    reeds_block.RETS = [std_RETS_val for i in range(len(SVALS))]

                if not isinstance(EIR, (Number, Iterable)) : #expand energy offsets to new s value ammount
                    EIR_vector = []
                    for x in reeds_block.EIR:
                        EIR_vector.append(x[0])

                    EIR_matrix = []
                    for z in EIR_vector:
                        EIR_matrix.append([z for i in range(len(SVALS))])

                    reeds_block.EIR = EIR_matrix

                reeds_block.RES = list(map(str, SVALS))
                reeds_block.NRES = len(SVALS)   #adjust number of Svals

            if (isinstance(RETS, (Iterable))):
                reeds_block.RETS = RETS

            if (isinstance(EIR, (Number, Iterable))):
                #set new EIR with 3 different types of input (single number, vector or matrix)
                EIR_matrix = []

                #single number
                if type(EIR) is int or type(EIR) is float:  # depends on SVALS and NRES
                    EIR_vector = [str(EIR) for x in range(int(reeds_block.NUMSTATES))]

                    for z in EIR_vector:
                        EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

                #vector or matrix
                elif type(EIR) is list and len(EIR) == int(reeds_block.NUMSTATES):
                    if len(EIR) == int(reeds_block.NUMSTATES) and type(EIR[0]) in [str, float, int]:
                        EIR_matrix = []
                        for z in EIR:
                            EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

                    elif(len(EIR) == int(reeds_block.NRES) and type(EIR[0]) is list and len(EIR[0]) == int(reeds_block.NUMSTATES)):
                        EIR_matrix = list(map(lambda x: list(map(str, x)), EIR))

                    else:
                        raise Exception(
                            "not enough EIR-vals for making REEDS Block. Got " + str(len(EIR)) + " for " + str(
                                reeds_block.NRES) + " SVals\n")
                reeds_block.EIR = EIR_matrix

            if (NRETRIAL ):
                reeds_block.NRETRIAL = NRETRIAL
            if (NREQUIL ):
                reeds_block.NREQUIL = NREQUIL
            if (CONT ):
                reeds_block.CONT = CONT

        elif(isinstance(reeds_block, blocks.REPLICA_EDS)):
            if(isinstance(REEDS, bool)):
                reeds_block.REEDS = REEDS

            if(isinstance(NUMSTATES, (Number, str))):
                reeds_block.NUMSTATES = NUMSTATES

            if(isinstance(SVALS, Iterable)): #edit SVALS
                SVALS = nice_s_vals(SVALS)

                reeds_block.RES = list(map(str, SVALS))
                reeds_block.NRES = len(SVALS)   #adjust number of Svals

                if(isinstance(EIR, (Number, Iterable))): #expand energy offsets to new s value ammount
                    # set new EIR with 3 different types of input (single number, vector or matrix)
                    EIR_matrix = []

                    # single number
                    if(isinstance(EIR, Number)):  # depends on SVALS and NRES
                        EIR_vector = [EIR for x in range(reeds_block.NUMSTATES)]

                        for z in EIR_vector:
                            EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

                    # vector or matrix
                    elif(isinstance(EIR, Iterable) and len(EIR) == int(reeds_block.NUMSTATES) and all([isinstance(x, Number) for x in EIR])):
                        for z in EIR:
                            EIR_matrix.append([z for i in range(int(reeds_block.NRES))])
                    else:
                        raise Exception(
                            "not enough EIR-vals for making REEDS Block. Got " + str(len(EIR)) + " for " + str(
                                reeds_block.NRES) + " SVals\n Number of states: "+str(reeds_block.NUMSTATES))
                    reeds_block.EIR = EIR_matrix
                else:
                    if(any([len(row)!= len(SVALS) for row in reeds_block.EIR])):
                        newEIR = []
                        for EIR_row in reeds_block.EIR:
                            newEIRrow = [EIR_row[0] for r in range(len(SVALS))]
                            newEIR.append(newEIRrow)
                        reeds_block.EIR = newEIR

            if(isinstance(EIR, (Number , Iterable))):
                EIR_matrix = []
                print(EIR)
                # single number
                if isinstance(EIR, Number):  # depends on SVALS and NRES
                    EIR_vector = [str(EIR) for x in range(reeds_block.NUMSTATES)]

                    for z in EIR_vector:
                        EIR_matrix.append([z for i in range(int(reeds_block.NRES))])

                # vector or matrix
                elif (isinstance(EIR, Iterable) and len(EIR) == int(reeds_block.NUMSTATES) and all(
                        [isinstance(x, (Number, str)) for x in EIR])):
                    for z in EIR:
                        EIR_matrix.append([float(z) for i in range(int(reeds_block.NRES))])
                else:
                    raise Exception(
                        "not enough EIR-vals for making REEDS Block. Got " + str(len(EIR)) + " for " + str(
                            reeds_block.NRES) + " SVals and states: "+str(reeds_block.NUMSTATES)+ "\n")
                reeds_block.EIR = EIR_matrix

            if (isinstance(NRETRIAL, (str, int))):
                reeds_block.NRETRIAL = NRETRIAL
            if (isinstance(NREQUIL, (str, int))):
                reeds_block.NREQUIL = NREQUIL
            if (isinstance(CONT, (str, int, bool))):
                reeds_block.CONT = CONT
            if (isinstance(EDS_STAT_OUT, (str, int, bool))):
                reeds_block.EDS_STAT_OUT = EDS_STAT_OUT

