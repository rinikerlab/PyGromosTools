"""
    Python implementation of the Gromos++ program ran_box which is used to generate randomized configurations for liquids (and gases)


    Author: Marc Lehner
"""
import warnings
import numpy as np
import copy
import random
import itertools as it


from pygromos.files.topology.top import Top
from pygromos.files.coord.cnf import Cnf

def ran_box(in_top_path:str, 
                    in_cnf_path:str, 
                    out_cnf_path:str= "", 
                    periodic_boundary_condition: str = "r", 
                    nmolecule:int = 1, 
                    dens:float = 1.0, 
                    threshold:float=None,
                    layer:bool = False, 
                    boxsize:float=None, 
                    fixfirst:bool = False, 
                    seed:float=None, 
                    _binary_name="ran_box", 
                    verbose=True, 
                    return_command_only=False)->str:

    top = Top(in_value=in_top_path)
    cnf = Cnf(in_value=in_cnf_path)
    cog = np.array(cnf.center_of_geometry())

    if(sum([len(cnf.residues[x]) for x in cnf.residues])>1):
        raise Exception("ran_box works only with one residue in the .cnf file!\nFound: "+str(cnf.get_residues()))

    #get volume and box length
    mol_mass = top.get_mass()
    volume = 1.66056 * nmolecule * mol_mass / dens
    box_length = volume**(1./3.)
    divider = int(np.ceil(nmolecule**(1./3.)))
    distance = box_length/float(divider)

    #calculate maxRandShift
    scale = 0.5 #scale can be manually decreased
    maxDist = 0
    for atom in copy.deepcopy(cnf).POSITION.content:
        pos = np.array([atom.xp, atom.yp, atom.zp])
        dis = np.linalg.norm(pos-cog)
        if dis > maxDist:
            maxDist = dis
    maxRandShift = (scale * distance) -maxDist
    if maxRandShift < 0:
        maxRandShift = 0
        if verbose:
            warnings.warn("Molecules might overlap! Check cnf manually or decrease the density")


    #create new cnf for return and set some attributes
    ret_cnf = copy.deepcopy(cnf)
    ret_cnf.POSITION.content = []
    if hasattr(ret_cnf, "LATTICESHIFTS"):
        delattr(ret_cnf, "LATTICESHIFTS")
    if hasattr(ret_cnf, "VELOCITY"):
        delattr(ret_cnf, "VELOCITY")
    if hasattr(ret_cnf, "STOCHINT"):
        delattr(ret_cnf, "STOCHINT")
    ret_cnf.GENBOX.pbc = 1
    ret_cnf.GENBOX.length = [box_length, box_length, box_length]
    ret_cnf.GENBOX.angles = [90,90,90]
    ret_cnf.TITLE.content = str(nmolecule) + " * " + cnf.POSITION.content[0].resName

    # add positions
    points = list(it.product(range(divider), range(divider), range(divider)))
    for ind, (xi, yi, zi) in enumerate(random.sample(points, nmolecule)):
        shift = np.array([xi*distance, yi*distance, zi*distance])
        cnf.rotate(alpha=random.uniform(0,360), beta=random.uniform(0,360), gamma=random.uniform(0,360))
        randomShift = np.array([random.uniform(-maxRandShift,maxRandShift),random.uniform(-maxRandShift,maxRandShift),random.uniform(-maxRandShift,maxRandShift)])

        for atom in copy.deepcopy(cnf).POSITION.content:
            pos = np.array([atom.xp, atom.yp, atom.zp])
            atom.xp, atom.yp, atom.zp = pos - cog + shift + randomShift
            atom.resID = ind+1
            atom.atomID += (ind * cnf.POSITION.content[-1].atomID)
            ret_cnf.POSITION.content.append(atom)
    
    ret_cnf.write(out_path=out_cnf_path)

    return out_cnf_path
