"""
File:            Class for trc files in pandas
Description:
    The pandas trajectory TRC class offers a easy method to process GROMOS's .trc files in python
    The trc files are parsed into an easy to use pandas dataframe

Author: Marc Thierry Lehner

TODO: add stride option to all member functions
TODO: add support for periodic boundary condition
TODO: add support for rdkit mol selector
TODO: add support for rdkit conformers

"""

#imports
import os
import pandas as pd
import numpy as np
from typing import TypeVar, Union
from pandas.core.base import DataError
import pygromos.files.trajectory._general_trajectory as traj
from pygromos.analysis import coordinate_analysis as ca

TrcType = TypeVar("Trc")
CnfType = TypeVar("Cnf")



class Trc(traj._General_Trajectory):
    _gromos_file_ending:str = "trc"
    def __init__(self, input_value: str or None, auto_save=True, stride:int=1, skip:int=0):
        super().__init__(input_value, auto_save=auto_save, stride=stride, skip=skip)


    def get_atom_pair_distance_series(self, atomI:int, atomJ:int, periodicBoundary=False) -> pd.DataFrame:
        """
        creates a series with the euclidian distance for the atom pair i-j for every time step

        Parameters
        ----------
        atomI : int
        atomJ : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
        """
        return self.database[["POS_"+str(atomI),"POS_"+str(atomJ)]].apply(lambda x: ca.calculate_distance(x[0], x[1]), axis=1)


    def get_atom_pair_distance_mean(self, atomI:int, atomJ:int) -> float:
        """
        creates the mean of the euclidian distances for every timestep between the atoms i and j

        Parameters
        ----------
        atomI : int
        atomJ : int

        Returns
        -------
        float
        """
        return self.get_atom_pair_distance_series(atomI=atomI, atomJ=atomJ).mean()


    def get_atom_movement_series(self, atomI:int, periodicBoundary=False) -> pd.DataFrame:
        """
        for every time step the difference to the last coordinates is calculated for an atom i
        Returns the euclidian difference to the last step as a numpy array

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
        """
        return (self.database["POS_"+str(atomI)]-self.database["POS_"+str(atomI)].shift())[1:]
    

    def get_atom_movement_length_series(self, atomI:int, periodicBoundary=False) -> pd.DataFrame:
        """
        for every timestep the total euclidian difference to the last time step is calculate for atom i

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
        """
        return self.get_atom_movement_series(atomI=atomI, periodicBoundary=periodicBoundary).apply(lambda x: np.sqrt(np.sum(x**2)))


    def get_atom_movement_length_mean(self, atomI:int, periodicBoundary=False) -> float:
        """
        The average euclidian movement between to consecutive timesteps is calculated for atom i

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        float
        """
        return self.get_atom_movement_length_series(atomI=atomI, periodicBoundary=periodicBoundary).mean()


    def get_atom_movement_length_total(self, atomI:int, periodicBoundary=False) -> float:
        """
        The total euclidian movement between to consecutive timesteps is calculated for atom i

        Parameters
        ----------
        atomI : int
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        float
        """
        return self.get_atom_movement_length_series(atomI=atomI, periodicBoundary=periodicBoundary).sum()


    def get_cog_movement_vector_series_for_atom_group(self, atoms:list, periodicBoundary=False) -> pd.DataFrame:
        """
        Calculate the movemnt vector between to consecutive time steps for the center of geometry defined by a list of atoms

        Parameters
        ----------
        atoms : list[int]
            the list of atoms defining the center of geometry
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
            numpy array with the movement vector for every time step
        """
        #create a local copy which can be modified
        temp_database = self.database.copy()
        # create the difference to the last time step for every atom i in the list atoms
        for i in atoms:
            temp_database["mix"+str(i)] = temp_database["POS_"+str(i)]-temp_database["POS_"+str(i)].shift()
        #calculate vector to cog
        return temp_database[["mix"+str(i) for i in atoms]][1:].sum(axis=1).apply(lambda x: x/len(atoms))

    def get_cog_movement_total_series_for_atom_group(self, atoms:list, periodicBoundary=False) -> pd.DataFrame:
        """
        Calculate the total movemnt between to consecutive time steps for the center of geometry defined by a list of atoms

        Parameters
        ----------
        atoms : list[int]
            the list of atoms defining the center of geometry
        periodicBoundary : bool, optional
            WIP

        Returns
        -------
        pd.DataFrame
            total movement (float) for every time step
        """
        return self.get_cog_movement_vector_series_for_atom_group(atoms=atoms, periodicBoundary=periodicBoundary).apply(lambda x: np.sqrt(np.sum(x**2)))

    def radial_distribution(self, atomsFrom, atomsTo) -> pd.DataFrame:
        pass


    def rmsd(self, ref_cnf: Union[int, TrcType]) -> pd.DataFrame:
        """Calculates the RootMeanSquareDeviation from a configuration (ref_cnf) to every frame in self

        Parameters
        ----------
        ref_cnf : trajectoy or frame (int) to use as reference
            This is the reference configuration. If a int (n) is provide the nth frame will be used

        Returns
        -------
        pd.DataFrame
            RMSD for every frame
        """
        if type(ref_cnf) == pd.DataFrame or type(ref_cnf) == pd.Series:
            if ref_cnf.ndim == 1:
                pos_mask=self.database.columns.str.startswith("POS_")
                to_sub = ref_cnf.iloc[pos_mask]
                if pos_mask.size != ref_cnf.size:
                    raise DataError("ref_cnf and Positons do not match")
        elif type(ref_cnf) == int:
            if ref_cnf <= self.database.ndim:
                pos_mask=self.database.columns.str.startswith("POS_")
                to_sub = self.database.iloc[ref_cnf, pos_mask]
            else:
                raise IndexError("ref_cnf value was recognized as integer but is too large")
        else:
            raise ValueError("ref_cnf type not supported")
        return self.database.iloc[:,pos_mask].sub(to_sub).apply(lambda x: ca.rms(x), axis=1)


    def get_pdb(self, cnf:str, exclude_resn=["SOLV"])->str:
        pdb_format = "ATOM  {:>5d}  {:<3}{:1}{:>3}  {:1}{:>3d}{:1}   {:>7.3f} {:>7.3f} {:>7.3f} {:>5}{:>6}{:<3}{:>2} {:>2d}"

        dummy_occupancy = dummy_bfactor = dummy_charge = 0.0
        dummy_alt_location = dummy_chain = dummy_insertion_code = dummy_segment = ""

        # 3) CONVERT FILE
        if( isinstance(self.TITLE, list)):
            TITLE = "\n".join(self.TITLE)
        else:
            TITLE = self.TITLE
        pdb_str = "TITLE " + TITLE + "\n"
        pos_cols = [col for col in self.database.columns if ("POS" in col)]
        for ind, time_step in self.database.iterrows():
            pos_lines = time_step[pos_cols]
            remark_line = "REMARK\t" + str(time_step["TIMESTEP_step"]) + "\t" + str(
                time_step["TIMESTEP_time"]) + "\nMODEL\n"
            frame_positions = []
            for ind, coord_set in enumerate(pos_lines):
                atom = cnf.POSITION[ind]
                if(atom.resName in exclude_resn):
                    continue
                frame_positions.append(
                    pdb_format.format(atom.atomID, atom.atomType, dummy_alt_location, atom.resName, dummy_chain, int(
                        atom.resID), dummy_insertion_code, coord_set[0] * 10, coord_set[1] * 10, coord_set[2] * 10,
                                      dummy_occupancy, dummy_bfactor,
                                      dummy_segment, atom.atomType, int(dummy_charge)))  # times *10 because pdb is in A
            pdb_str+= remark_line + "\n".join(frame_positions) + '\nENDMDL\n'
        return pdb_str


    def write_pdb(self, out_path:str, cnf_file:str):
            """
                This function converts the atom POS db of the traj into a pdb traj.

            Parameters
            ----------
            atoms : t.List[Atom]
                List of atoms
            Returns
            -------
            t.List[str]
                 pdb strings of that molecule
            """
            # 1) INPUT PARSING
            from pygromos.files.coord.cnf import Cnf # avoid circular import
            if(isinstance(cnf_file, str)):
                cnf = Cnf(cnf_file)
            elif(isinstance(cnf_file, Cnf)):
                cnf = cnf_file
            else:
                raise ValueError("Did not understand the Value of cnf_file. Must be str to a cnf file or a cnf_file.")

            if(isinstance(out_path, str)):
                if(os.path.exists(os.path.dirname(out_path))):
                    out_file = open(out_path, "w")
                else:
                    raise IOError("Could not find directory to write to: "+str(os.path.dirname(out_path)))
            else:
                raise ValueError("Did not understand the Value of out_path. Must be str.")

            # 2) CONSTUCT PDB BLOCKS
            # ref: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
            pdb_format = "ATOM  {:>5d}  {:<2}{:1}{:>3}  {:1}{:>3d}{:1}   {:>7.3f}{:>7.3f}{:>7.3f}{:>5}{:>6}{:<3}{:>2} {:>2d}"
            dummy_occupancy = dummy_bfactor = dummy_charge = 0.0
            dummy_alt_location = dummy_chain = dummy_insertion_code = dummy_segment = ""

            # 3) CONVERT FILE
            # TODO: Inefficient!
            out_file.write( "TITLE " + self.TITLE + "\n")
            pos_cols = [col for col in self.database.columns if("POS" in col)]
            for ind, time_step in self.database.iterrows():
                pos_lines = time_step[pos_cols]
                remark_line="REMARK\t"+str(time_step["TIMESTEP_step"])+"\t"+str(time_step["TIMESTEP_time"])+"\nMODEL\n"
                frame_positions = []
                for ind ,coord_set in enumerate(pos_lines):
                    atom =cnf.POSITION[ind]
                    frame_positions.append(
                        pdb_format.format(atom.atomID, atom.atomType, dummy_alt_location, atom.resName, dummy_chain, int(
                            atom.resID), dummy_insertion_code, coord_set[0]*10, coord_set[1]*10, coord_set[2]*10, dummy_occupancy, dummy_bfactor,
                        dummy_segment, atom.atomType, int(dummy_charge))) #times *10 because pdb is in A
                out_file.write(remark_line+"\n".join(frame_positions)+'\nENDMDL\n')
            out_file.write('\nEND')
            out_file.close()
            return out_path


    def visualize(self, cnf:CnfType):
        from pygromos.visualization.coordinates_visualization import show_coordinate_traj
        return show_coordinate_traj(self, cnf=cnf)


    def cog_reframe(self, cnf:CnfType, index_list:list=[1]):
        # create mask for cog calculation
        col_list = [x for x in self.database.columns if ("POS" in x)]

        # calculate average box size
        grid = self.database["length"].mean()
        grid = np.array( cnf.GENBOX.length) / 2

        #print("Grid", grid)

        # cog calculation: select POS -> apply pbc -> average all positions
        pbc_pos = self.database[col_list].applymap(lambda x: ca._periodic_distance(x, grid))
        cog = pbc_pos.sum(axis=1) / len(col_list)
        #print("COG:", cog)

        # box center
        boxCenter = self.database["length"] / 2
        #print("box_center:", boxCenter.shape, boxCenter[0])

        # shift all positions
        posList = [x for x in self.database.columns if x.startswith('POS_')]
        for ind, idx in enumerate(posList):
            self.database[idx] = self.database[idx].apply(lambda x: self._periodic_distance(x, grid))  # + boxCenter

"""
    @classmethod
    def cnfs_to_trc(cls, cnfs:Iterable[str], title=None, start_time:float=0.0, dt:float=0.002 ):
        return cls(input=cls._cnf_to_trc_dict(cnfs, start_time=start_time, dt=dt), title=title)

    @classmethod
    def _cnf_to_trc_dict(cls, cnfs:Iterable[str], start_time:float=0.0, dt:float=0.002)->List[Dict]:

        time = start_time
        out_frames_list:List[Dict] = []
        for steps, cnf_path in enumerate(cnfs):
            cnf = Cnf(cnf_path)
            frame_traj = {'steps': steps, 'time': time, "POSITIONRED": []}
            block_str = ""
            for atom in sorted(cnf.POSITION, key=lambda at: at.atomID):
                if (atom.atomID % 10 == 0): block_str += "#\t" + str(atom.atomID) + "\n"
                block_str += "\t" + str(atom.xp) + "\t" + str(atom.yp) + "\t" + str(atom.zp) + "\n"
            frame_traj.update({"POSITIONRED": block_str})
            frame_traj.update({"GENBOX": str(cnf.GENBOX).replace("GENBOX\n", "").replace("END\n", "")}) #a bit hacky

            out_frames_list.append(frame_traj)
            time += dt

        return out_frames_list
"""
