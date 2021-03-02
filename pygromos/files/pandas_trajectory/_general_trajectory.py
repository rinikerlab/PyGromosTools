"""
FUNCTIONLIB:            trajectory template with pandas
Description:
    in this parent class for all trajectories the general parsing to a pandas database is handeled.
    This class is intended as quick data evaluation in python without gromos++ programs like ene_ana

    Specific block structures are handeled in the respectve classes (trc/tre/trf/...)

Author: Marc Thierry Lehner

Remark: WIP

###########################
ATTENTION: NOT FULLY TESTET
ATTENTION: WORK IN PROGRESS
###########################

TODO: MTL: implement, check, test

"""

#imports
import warnings, glob, importlib, collections, os
import numpy
import pandas

import pygromos.files.pandas_trajectory.trajectory_blocks as blocks

class _General_Trajectory():
    def __init__(self, input:(str or None)):
        if input == None:
            raise "Empty construction not impokemented yet"
        elif type(input) == str:
            self.database = self.read_trajectory(input=input)
        else:
            raise "Constructor not found"
    
    def _raw_read_trajectory(self, input:str) -> collections.defaultdict:
        #define temp storage
        header = {}
        data = []
        dataInTimestep = {}
        block = []
        blockname = ''
        
        #set contorl bool
        isInTimeStep = False
        isInBlock = False

        # start to read the trajectory
        with open(input, 'r') as infile:
            for line in infile:
                if isInTimeStep:
                    if isInBlock:
                        if not line.strip().startswith("END"):
                            block.append(line)
                        else:
                            isInBlock = False
                            dataInTimestep.update({blockname:block})                        
                    else:
                        if line.strip().startswith('#') or line.strip() == '':
                            continue
                        blockname = line.strip()
                        block = []
                        isInBlock = True
                        if blockname.startswith("TIMESTEP"):
                            data.append(dataInTimestep)
                            dataInTimestep = {}
                else:
                    if isInBlock:
                        if not line.strip().startswith("END"):
                            block.append(line)
                        else:
                            isInBlock = False
                            header.update({blockname:block})
                    else:
                        if line.strip().startswith('#') or line.strip() == '':
                            continue
                        blockname = line.strip()
                        block = []
                        isInBlock = True
                        if blockname.startswith("TIMESTEP"):
                            isInTimeStep = True
            if isInTimeStep:
                #final time step finish since no smarter way to detect end of Timestep
                data.append(dataInTimestep)
            else:
                raise "No timestep found"
        return {"header":header, "body":data}

    def read_trajectory(self, input:str) -> pandas.DataFrame:
        if (not os.path.exists(input)):
            raise IOError("Could not find File: ", input)
        else:
            table = []
            data = self._raw_read_trajectory(input=input)
            header = data["header"]
            body = data["body"]
            for key, block in header.items():
                setattr(self, key, block)
            for time_step_entry in body:
                table_entry = {}
                for blocktitle, block in time_step_entry.items():
                    if not hasattr(blocks, blocktitle):
                        raise "No trajectory block found named: "+blocktitle
                    tmp_block = getattr(blocks, blocktitle)(block)
                    table_entry.update(tmp_block.to_dict())
                table.append(table_entry)
        return(pandas.DataFrame.from_dict(table))


        
        
                    

                
                




