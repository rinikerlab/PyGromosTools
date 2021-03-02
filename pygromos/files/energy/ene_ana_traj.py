"""
    NOT a doc string
"""
import os
import warnings
from typing import List, Dict, Any

from pygromos.utils import bash as bash
import pandas as pd

class ene_ana_property_traj():
    """
        warning: CONCIDER DROPPing this class
    """
    name:str
    in_path:List[str] or None
    path:str or List[str]
    properties:Dict[str, List[Any]] #any property/ies
    time:List[float] #in ps

    def __init__(self, in_path:str=None, name:str=None, properties:Dict[str, List[Any]]=None, time:List[float]=None, verbose:bool=False)->None:
        if(in_path!= None):
            self.name, self.properties,  self.time = self.parse_property_traj(in_path, verbose=verbose)
            self.in_path = [in_path]

        elif(all([x!=None for x in [name, properties, time]])):
            self.name=name
            self.properties=properties
            self.time=time
            self.in_path = None
        else:
            raise IOError("could not find correctd constructor parameters!")

    def __str__(self)->str:
        msg = "# " + self.name + "\n"
        msg+="# time\t" + "\t".join([protperty for protperty in sorted(self.properties)]) + "\n"
        for ind, time in enumerate(self.time):
            msg+=" " + str(time) + "\t" + "\t".join([str(self.properties[property][ind]) for property in sorted(self.properties)])+"\n"
            pass
        return msg

    @staticmethod
    def parse_property_traj(path:str, verbose:bool=False)->(str, Dict[str, List[Any]], List[float]):
        e_file = open(path, "r")  # open file
        e_lines = e_file.readlines()
        e_comments = list(filter(lambda x: x.startswith("#"), e_lines))
        e_header = e_comments[-1].strip().split()[1:]  # get headers
        e_lines = filter(lambda x: not x.startswith("#"), e_lines)  # filter comments
        e_fields = map(lambda x: x.strip().split(), e_lines)  # split lines and convert to float
        e_fields = list(map(lambda x: tuple(map(float, x)), e_fields))  # split lines and convert to float

        if verbose: print("FILE: Header: ", e_header)

        #name
        if(len(e_comments)>1):
            name = e_comments[0].strip()
        else:
            name="unknown"

        #TIME
        time = list(map(lambda x: x[0], e_fields))  # extract time
        first = False

        # collect potentials
        num_properties = len(e_header[1:])
        properties = {}
        property_values = [list(map(lambda x: x[i], e_fields)) for i in range(1, num_properties + 1)]
        if (all([len(property_values[i]) == len(time) for i in range(num_properties)])):
            properties.update({e_header[i + 1]: property_values[i] for i in range(num_properties)})
        else:   #CHECK if same length as time.
            raise ValueError("The time length of state 1 is not fitting to the ammount of Energy values in traj name: " + str(name) +
                             "\n\t " + str(len(time)) + " vs. " + str(len(property_values)))

        return name, properties, time

    def write(self, out_path:str, verbose:bool=False, safe:bool=True, remove_input_files:bool=False)->(None or str):
        # write Output_file
        if (not os.path.isdir(os.path.dirname(out_path))):
            raise ValueError("Could not find output folder: " + out_path + " \n Could not write out energy tra: \n" + out_path)

        if(os.path.exists(out_path) and safe):
            warnings.warn("Found File already! Aborting write out.\n Path: "+out_path)
            return None

        if verbose: print("Writing out: "+out_path)

        #write out!
        out_file = open(out_path, "w")
        if verbose: print("# " + self.name + "\n")
        out_file.write("# " + self.name + "\n")  # print title
        if verbose: print("# time\t" + "\t".join([protperty for protperty in sorted(self.properties)]) + "\n")
        out_file.write("# time\t" + "\t".join([protperty for protperty in sorted(self.properties)]) + "\n")  # make header

        ##write_data:
        for ind, time in enumerate(self.time):
            if(verbose): print(" " + str(time) + "\t" + "\t".join(
                [str(self.properties[property][ind]) for property in sorted(self.properties)]))
            out_file.write(" " + str(time) + "\t" + "\t".join([str(self.properties[property][ind]) for property in sorted(self.properties)]) + "\n")
        out_file.close()
        self.path=out_path

        #remove unconcatenated input_files
        if(remove_input_files):
            self.remove_input_files(verbose=verbose)
        return out_path

    def remove_input_files(self, verbose:bool=False):
        if(type(self.in_path) == str):
            if (self.path == self.in_path):            #safety checkup
                warnings.warn("Did not remove in_dir as it is also current out_path!")
                self.in_path = None
                return
            if verbose: print("DELETE FILE: " + self.in_path)
            bash.remove_file(self.in_path)

        elif(type(self.in_path) == list):
            if (self.path in self.in_path):            #safety checkup
                warnings.warn("Found current out_path in in_dir List, removed it from there, continue cleaning.")
                pos = self.in_path.index(self.path)
                del self.in_path[pos]

            for in_p in self.in_path:
                if(type(in_p) == str):
                    if verbose: print("DELETE FILE: " + in_p)
                    bash.remove_file(in_p)
                else:
                    raise ValueError("in_p is not a str but: "+str(type(in_p))+" with value "+str(in_p))
        elif(type(self.in_path) == None):
            warnings.warn("No in_paths found!")
        else:
            raise Exception("Could not understand in_dir attribute type, during cleanup!\n in_paths was: "+str(self.in_path))
        self.in_path = None

    def add_property_file(self, in_path:str=None, in_traj_obj:object=None, choose_properties:List[str]=None,
                          overwrite_time:bool=False, overwrite_property:bool=False, verbose:bool=False)->None:
        #parse
        if(in_path):
            _, properties, time = self.parse_property_traj(in_path, verbose=verbose)
            if(self.in_path == None):
                self.in_path == [in_path]
            else:
                self.in_path.append(in_path)

        elif(in_traj_obj):
            properties, time = (in_traj_obj.properties, in_traj_obj.time)
        else:
            raise IOError("For add_properties give either in_obj or in_dir!")

        #chosing subset of potentials???
        if(choose_properties != None):

            if(any([not prop in properties for prop in  choose_properties])):
                raise IOError("Could not find chosen property arg in found potentials!\n"
                              " missing chosenProp:\n"+"\t".join([prop for prop in  choose_properties if(not prop in properties)]+"")+"\n"
                              " found props in file: \n"+"\t".join(list(properties.keys())))
            if(verbose):print("filter for selected potentials: "+"\t".join(properties.keys())+"\n"
                              "------------------------------> "+"\t".join(choose_properties))
            chosen_properties = {key: properties[key] for key in choose_properties}
        else:
            chosen_properties=properties

        #check different props.
        if(len(time)==len(self.time)):
            if(all([t1==t2 for t1, t2 in zip(self.time, time)]) or overwrite_time):
                if(any([prop in self.properties for prop in chosen_properties]) and not overwrite_property):
                    raise IOError("Found in both files the same Property! \n "
                                  "chosen Props: \n"+"\t".join([prop for prop in chosen_properties if(prop in properties)])+"\n"
                                  "already existing potentials: \n"+"\t".join(list(self.properties.keys())))
                else:
                    if(verbose):print("Adding potentials: "+"\t".join(chosen_properties.keys()))
                    [self.properties.update({key:value}) for key, value in chosen_properties.items()]
            else:
                raise IOError("TIME-Values are not the same in the files! \nuncorrect tuples: \n"+"\n".join([str(t1)+"\t"+str(t2) for t1, t2 in zip(self.time, time) if(t1!=t2)]))
        else:
            raise IOError("There are not the same ammount of data points in both files! (you could implement a mapping) :)\n"
                          "obj-time: "+str(len(self.time))+"\n"
                          "file-time: "+str(len(time))+"\n")

    def add_property_files(self, in_trajs:List[str or object], choose_properties:List[str]=None,
                           overwrite_time:bool=False, overwrite_property:bool=False, verbose:bool=False)->None:
        for in_traj in in_trajs:
            if(type(in_traj) == str):
                self.add_property_file(in_path=in_traj, choose_properties=choose_properties,
                                       overwrite_time=overwrite_time, overwrite_property=overwrite_property,verbose=verbose)
            else:
                self.add_property_file(in_traj_obj=in_traj, choose_properties=choose_properties,
                                       overwrite_time=overwrite_time, overwrite_property=overwrite_property, verbose=verbose)
        pass

    def reduce_data(self, num_data_points:int, verbose:bool=False):

        if(num_data_points>len(self.time)):
            raise ValueError("The number, which the data should be reduced to, is bigger than present data.\n num_points" + str(num_data_points) + "\tpresent: " + str(len(self.time)))
        elif(len(self.time)%num_data_points!=0):
            raise ValueError("The number, which the data should be reduced to, is not cleanly dividable by present data.\n num_data_points%present_data!= 0 but " + str(num_data_points % len(self.time)))

        dstep = int(len(self.time)/num_data_points)
        if(verbose): print("take each "+str(dstep)+" step")

        self.time = self.time[::dstep]

        for prop in self.properties:
            if(verbose): print("REDUCE "+str(prop))
            tmp = self.properties[prop][::dstep]
            self.properties.update({prop:tmp})



class new_ene_ana_property_traj():
    """

    """
    name:str
    in_path:List[str] or None
    path:str or List[str]
    properties:Dict[str, List[Any]] #any property/ies
    time:List[float] #in ps

    def __init__(self, in_path:str=None, name:str=None,
                 sep:str="\t", header_lines:int=1, verbose:bool=False)->None:
        if(isinstance(in_path, str)):
            self.data = pd.read_csv(in_path, header=header_lines, sep=sep)
            self.in_path = [in_path]

        else:
            super().__init__()

        if(isinstance(name, str)):
            self.name = name
            self.in_path = []

        #else:
        #    raise IOError("could not find correctd constructor parameters!")

    @staticmethod
    def parse_property_traj(path:str, verbose:bool=False)->(str, Dict[str, List[Any]], List[float]):
        e_file = open(path, "r")  # open file
        e_lines = e_file.readlines()
        e_comments = list(filter(lambda x: x.startswith("#"), e_lines))
        e_header = e_comments[-1].strip().split()[1:]  # get headers
        e_lines = filter(lambda x: not x.startswith("#"), e_lines)  # filter comments
        e_fields = map(lambda x: x.strip().split(), e_lines)  # split lines and convert to float
        e_fields = list(map(lambda x: tuple(map(float, x)), e_fields))  # split lines and convert to float

        if verbose: print("FILE: Header: ", e_header)

        #name
        if(len(e_comments)>1):
            name = e_comments[0].strip()
        else:
            name="unknown"

        #TIME
        time = list(map(lambda x: x[0], e_fields))  # extract time
        first = False

        # collect potentials
        num_properties = len(e_header[1:])
        properties = {}
        property_values = [list(map(lambda x: x[i], e_fields)) for i in range(1, num_properties + 1)]
        if (all([len(property_values[i]) == len(time) for i in range(num_properties)])):
            properties.update({e_header[i + 1]: property_values[i] for i in range(num_properties)})
        else:   #CHECK if same length as time.
            raise ValueError("The time length of state 1 is not fitting to the ammount of Energy values in traj name: " + str(name) +
                             "\n\t " + str(len(time)) + " vs. " + str(len(property_values)))

        return name, properties, time

    def write(self, out_path:str, remove_input_files:bool=False, safe:bool=False)->str:
        if(safe and os.path.exists(out_path)):
            raise Exception("File already exists: "+ out_path)
        else:
            self.data.to_csv(out_path, sep="\t")

            if(remove_input_files):
                self.remove_input_files()
        return out_path

    def remove_input_files(self, verbose:bool=False):
        if(type(self.in_path) == str):
            if (self.path == self.in_path):            #safety checkup
                warnings.warn("Did not remove in_dir as it is also current out_path!")
                self.in_path = None
                return
            if verbose: print("DELETE FILE: " + self.in_path)
            bash.remove_file(self.in_path)

        elif(type(self.in_path) == list):
            if (self.path in self.in_path):            #safety checkup
                warnings.warn("Found current out_path in in_dir List, removed it from there, continue cleaning.")
                pos = self.in_path.index(self.path)
                del self.in_path[pos]

            for in_p in self.in_path:
                if(type(in_p) == str):
                    if verbose: print("DELETE FILE: " + in_p)
                    bash.remove_file(in_p)
                else:
                    raise ValueError("in_p is not a str but: "+str(type(in_p))+" with value "+str(in_p))
        elif(type(self.in_path) == None):
            warnings.warn("No in_paths found!")
        else:
            raise Exception("Could not understand in_dir attribute type, during cleanup!\n in_paths was: "+str(self.in_path))
        self.in_path = None

    def add_property_file(self, in_path:str=None, in_traj_obj:object=None, choose_properties:List[str]=None,
                          overwrite_time:bool=False, overwrite_property:bool=False, verbose:bool=False)->None:
        #parse
        if(in_path):
            _, properties, time = self.parse_property_traj(in_path, verbose=verbose)
            if(self.in_path == None):
                self.in_path == [in_path]
            else:
                self.in_path.append(in_path)

        elif(in_traj_obj):
            properties, time = (in_traj_obj.properties, in_traj_obj.time)
        else:
            raise IOError("For add_properties give either in_obj or in_dir!")

        #chosing subset of potentials???
        if(choose_properties != None):

            if(any([not prop in properties for prop in  choose_properties])):
                raise IOError("Could not find chosen property arg in found potentials!\n"
                              " missing chosenProp:\n"+"\t".join([prop for prop in  choose_properties if(not prop in properties)]+"")+"\n"
                              " found props in file: \n"+"\t".join(list(properties.keys())))
            if(verbose):print("filter for selected potentials: "+"\t".join(properties.keys())+"\n"
                              "------------------------------> "+"\t".join(choose_properties))
            chosen_properties = {key: properties[key] for key in choose_properties}
        else:
            chosen_properties=properties

        #check different props.
        if(len(time)==len(self.time)):
            if(all([t1==t2 for t1, t2 in zip(self.time, time)]) or overwrite_time):
                if(any([prop in self.properties for prop in chosen_properties]) and not overwrite_property):
                    raise IOError("Found in both files the same Property! \n "
                                  "chosen Props: \n"+"\t".join([prop for prop in chosen_properties if(prop in properties)])+"\n"
                                  "already existing potentials: \n"+"\t".join(list(self.properties.keys())))
                else:
                    if(verbose):print("Adding potentials: "+"\t".join(chosen_properties.keys()))
                    [self.properties.update({key:value}) for key, value in chosen_properties.items()]
            else:
                raise IOError("TIME-Values are not the same in the files! \nuncorrect tuples: \n"+"\n".join([str(t1)+"\t"+str(t2) for t1, t2 in zip(self.time, time) if(t1!=t2)]))
        else:
            raise IOError("There are not the same ammount of data points in both files! (you could implement a mapping) :)\n"
                          "obj-time: "+str(len(self.time))+"\n"
                          "file-time: "+str(len(time))+"\n")

    def add_property_files(self, in_trajs:List[str or object], choose_properties:List[str]=None,
                           overwrite_time:bool=False, overwrite_property:bool=False, verbose:bool=False)->None:
        for in_traj in in_trajs:
            if(type(in_traj) == str):
                self.add_property_file(in_path=in_traj, choose_properties=choose_properties,
                                       overwrite_time=overwrite_time, overwrite_property=overwrite_property,verbose=verbose)
            else:
                self.add_property_file(in_traj_obj=in_traj, choose_properties=choose_properties,
                                       overwrite_time=overwrite_time, overwrite_property=overwrite_property, verbose=verbose)
        pass

    def reduce_data(self, num_data_points:int, verbose:bool=False):

        if(num_data_points>len(self.time)):
            raise ValueError("The number, which the data should be reduced to, is bigger than present data.\n num_points" + str(num_data_points) + "\tpresent: " + str(len(self.time)))
        elif(len(self.time)%num_data_points!=0):
            raise ValueError("The number, which the data should be reduced to, is not cleanly dividable by present data.\n num_data_points%present_data!= 0 but " + str(num_data_points % len(self.time)))

        dstep = int(len(self.time)/num_data_points)
        if(verbose): print("take each "+str(dstep)+" step")

        self.time = self.time[::dstep]

        for prop in self.properties:
            if(verbose): print("REDUCE "+str(prop))
            tmp = self.properties[prop][::dstep]
            self.properties.update({prop:tmp})
