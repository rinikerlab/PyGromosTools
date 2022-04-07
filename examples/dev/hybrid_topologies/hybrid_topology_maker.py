import os, copy, glob

from datetime import datetime

import numpy as np
import rdkit

import pygromos
from pygromos.files.blocks import topology_blocks as blocks
from pygromos.files.topology.top import Top
from pygromos.files.coord.cnf import Cnf
from pygromos.files.blocks.coord_blocks import atomP


def reduceConformations(path_cnf, out_path, contains_protein=False):
    """
    This function will make a cnf for each possible ligand (alone in water or in complex)
    based on the given re-eds dual topology type cnf. 
    
    These reduced conformations can be used for single ligand simulations, as well as 
    to reconstruct the Hybrid conformation. 
    
    Parameters
    ----------
        path_cnf: str
            path to input cnf
        out_path: str
            path to output file
        contains_protein: bool
    
    Returns
    --------
        None
    
    """
    
    # Open the cnf containing everything. 
    cnf = Cnf(path_cnf)

    # Count number of ligands
    residues = cnf.get_residues() 
    if 'WAT' in residues: del residues['WAT']
    if 'SOLV' in residues: del residues['SOLV']
    num_ligs = 0
    for key, value in residues.items():
        if 'LI' in key: num_ligs += 1
    
    # Make a cnf for each single ligand solvated 

    for i in range(num_ligs):
        # 1: Make a copy
        tmp_cnf = copy.deepcopy(cnf)

        keep_res_at_id1 = False

        # 2: Remove all other ligands
        for j in range(num_ligs):
            if i == j: 
                keep_res_at_id1 = True
                continue
            
            # delete_residue renumbers residues so we need to keep track of which one to remove
            if (keep_res_at_id1):
                tmp_cnf.delete_residue(resID=(2))
            else:
                tmp_cnf.delete_residue(resID=(1))
            
        # 3: Print results to a file
        if contains_protein:
            tmp_cnf.write(out_path = out_path +'/complex_' + str(i+1) + '.cnf')
            tmp_cnf.write_pdb(out_path = out_path +'/complex_' + str(i+1) + '.pdb')
        else:
            tmp_cnf.write(out_path = out_path +'/ligand_' + str(i+1) + '.cnf')
            tmp_cnf.write_pdb(out_path = out_path +'/ligand_' + str(i+1) + '.pdb')

    return None


#  We will keep track of all of the data of PerturbedAtoms
#  in this class, by making a list of PerturbedAtoms (called atom_mappings)
#  which will then allow us to construct the cnf and ptp files
class PerturbedAtom:
    
    def __init__(self, atom:str , init_lig:int, init_id:int, new_id:int):
        """
        Create object with data given as input
        """
        self.atom = atom
        self.init_lig = init_lig
        self.init_id = init_id
        self.new_id = new_id

    def __str__(self):
        return ('atom: ' + self.atom + '   ' + 'initial ligand: ' + str(self.init_lig)
              + '  initial id: ' + str(self.init_id) + '   new_id: ' + str(self.new_id))
    
def find_new_id(atom_mapping, init_lig, init_id):
        """
            This function will return the new ID of the atom which belong 
            to the initial ligand init_lig and and an initial id init_id, 
            from the list of atom mappings 
            
            Parameters
            ----------
                mappings: List [PerturbedAtom]
                    
        """
        
        for am in atom_mapping:
            if am.init_lig == init_lig and am.init_id == init_id: return am.new_id
        
        # if we reach this we have a problem
        raise Exception("We did not find the atom you were looking for in the list\n" +
                        'init_lig: ' + str(init_lig) + ' and init_id: ' + str(init_id))

#
#
#
def findAtomsInCoreAndRGroup(topo, connection_point): 
    """
        This function will find the list of atoms in a molecule
        which belong to the "core" region, and those which belong 
        to the "-R group"
        
        This function even works when there are multiple connecting points,  
        (i.e. if multiple -R groups exist)        
 
        Parameters
        ----------
        topo: pygromos topology 
        
        connection_point: (int, int) or (List[int], List[int])
            first element is the atom (or list of atoms) in the core region 
            second element is the atom (or list of atoms) in the rgroup
        Returns
        ---------
            (core_atoms, rgroup_atoms): (List[int], List[int])
            atom ids of the core region and rgroup in the topology
    """
    
    # 1: Make the list of all bonds in the molecule
    
    bonds = [] # will be a list of tuples
    
    for bond in topo.BOND.content:
        bonds.append( (bond.IB, bond.JB) )
    for bond in topo.BONDH.content:
        bonds.append( (bond.IB, bond.JB) )
    
    num_atoms = topo.SOLUTEATOM.content[-1].ATNM
    
    # 2: Initialize results data structure
    
    if isinstance(connection_point[0], int): 
        core_atoms = [connection_point[0]]
    else:
        core_atoms = connection_point[0]
    
    if isinstance(connection_point[1], int): 
        rgroup_atoms = [connection_point[1]]
    else:
        rgroup_atoms = connection_point[1]
    
    # 3: Iterate over the list until all atoms have been categorized
    
    niter = 0
    
    while len(core_atoms) + len(rgroup_atoms) != num_atoms:
        for bond in bonds:            
            # If its the connecting bond, this needs to be skipped
            if bond[0] in core_atoms and bond[1] in rgroup_atoms or \
               bond[1] in core_atoms and bond[0] in rgroup_atoms: continue  
            
            if   bond[0] in core_atoms and bond[1] not in core_atoms: core_atoms.append(bond[1])
            elif bond[1] in core_atoms and bond[0] not in core_atoms: core_atoms.append(bond[0])
            
            if   bond[0] in rgroup_atoms and bond[1] not in rgroup_atoms: rgroup_atoms.append(bond[1])
            elif bond[1] in rgroup_atoms and bond[0] not in rgroup_atoms: rgroup_atoms.append(bond[0])
        
        niter += 1
        if niter > 500: break # This is just to kill the unfinite loop
        
    return (sorted(core_atoms), sorted(rgroup_atoms))


# This file contains all of the functions used to build hybrid topologies
# the functions were moved here to mke the jupyter notebook lighter.

def addHybridAtoms(core_top, new_lig_top, new_lig_rgroup, atom_mappings):
    """
        This function will add all of the atoms from a new 
        ligand into the topology with the core. 
    
        Parameters
        ----------
            core_top: topology
                Core gromos topology object which will be updated (atoms added to this object)
            new_lig_top: topology
                Gromos topology object of the other ligand for which we will add atoms to core_top
            new_lig_rgroup: List [int]
                List of atomIDs in new_lig_top to add to core_top
            atom_mappings: List [HybridAtoms]
                List of HybridAtoms object keeping track of the atom renumbering
            
        Returns
        --------
            core_top, atom_mappings: Topology, List[HybridAtoms]
                core_top is the updated topology object
                atom_mappings is the list keeping track of the atom renumbering
    """
    atnmShift = core_top.SOLUTEATOM.content[-1].ATNM #Number of atoms found in main top. Shift secondary top atoms accordingly
    mresShift = core_top.SOLUTEATOM.content[-1].MRES #Number of molecules found in main top.

    # Loop over all atoms in the new topology and add if it is in new_lig_rgroup
    atom_counter = 0
    
    for atom in new_lig_top.SOLUTEATOM.content:
        if atom.ATNM not in new_lig_rgroup: continue
        atom_counter += 1
        core_top.add_new_soluteatom(ATNM = atnmShift + atom_counter,
                                    MRES = mresShift + atom.MRES,
                                    PANM = atom.PANM,
                                    IAC = atom.IAC,
                                    MASS = atom.MASS,
                                    CG = atom.CG, 
                                    CGC = 0, # we will make the entire perturbed region one charge group.
                                    INE = [], # these will be filled in after all ligands have been added!
                                    INE14 = [])
        
        # Keep track of the atoms that were added in the mapping list. 
        atom_mappings.append(PerturbedAtom(atom = atom.PANM, 
                                           init_lig = mresShift + atom.MRES, 
                                           init_id  = atom.ATNM,
                                           new_id   = atnmShift + atom_counter))
    return core_top, atom_mappings


def addHybridBonds(core_top, new_lig_top, new_lig_rgroup, atom_mappings):
    """
        This function will add all of the bonds from a new 
        ligand into the topology with the core. 
        
        This includes the bonds that the -R group forms with itself
        + the one connecting bond between the core region and the R group.
    
        Parameters
        ----------
            core_top: topology
                Core gromos topology object which will be updated (atoms added to this object)
            new_lig_top: topology
                Gromos topology object of the other ligand for which we will add atoms to core_top
            new_lig_rgroup: List [int]
                List of atomIDs in new_lig_top to add to core_top
            atom_mappings: List [HybridAtoms]
                List of HybridAtoms object keeping track of the atom renumbering
            
        Returns
        --------
            core_top, atom_mapping: Topology, List[HybridAtoms]
                core_top is the updated topology object
                atom_mappings is the list keeping track of the atom renumbering
    """
    
    # Since this will be called just after adding the bonds, 
    # we assume that the number of residues in the topology equals the residue we are currently working for.
    # NOTE: might have to change this if I make the core different than the -R groups? 
    mres = core_top.SOLUTEATOM.content[-1].MRES
    # Debugging statement to see things more clearly in the file. 
    
    debugging_mode = False
    
    if debugging_mode:
        newBond = blocks.top_bond_type(IB=99999, JB=99999, ICB=99999)
        core_top.BOND.content.append(newBond)
        core_top.BOND.NBON += 1
    
    for bond in new_lig_top.BOND.content:
        if (bond.IB not in new_lig_rgroup) and (bond.JB not in new_lig_rgroup): continue
        
        # note: we need to atom_mappings to be properly set here.
        newBond = blocks.top_bond_type(IB=find_new_id(atom_mappings, mres, bond.IB),
                                       JB=find_new_id(atom_mappings, mres, bond.JB),
                                       ICB=bond.ICB)    
        # Append this new bond to the topology
        core_top.BOND.content.append(newBond)
        core_top.BOND.NBON += 1
        
    # Do exactly the same thing for bonds containing Hydrogens
    
    if debugging_mode:
        newBond = blocks.top_bond_type(IB=99999, JB=99999, ICB=99999)
        core_top.BONDH.content.append(newBond)
        core_top.BONDH.NBONH += 1
    
    for bondh in new_lig_top.BONDH.content:
        if (bondh.IB not in new_lig_rgroup) and (bondh.JB not in new_lig_rgroup): continue
        
        newBond = blocks.top_bond_type(IB=find_new_id(atom_mappings, mres, bondh.IB),
                                       JB=find_new_id(atom_mappings, mres, bondh.JB),
                                       ICB=bond.ICB)    
        # Append this new bond to the topology
        core_top.BONDH.content.append(newBond)
        core_top.BONDH.NBONH += 1
    
    return core_top, atom_mappings

def addHybridAngles(core_top, new_lig_top, new_lig_rgroup, atom_mappings):
    """
        This function will add all of the angles from a new 
        ligand into the topology with the core. 
        
        This includes the bonds that the -R group forms with itself
        + the one connecting bond between the core region and the R group.
    
        Parameters
        ----------
            core_top: topology
                Core gromos topology object which will be updated (atoms added to this object)
            new_lig_top: topology
                Gromos topology object of the other ligand for which we will add atoms to core_top
            new_lig_rgroup: List [int]
                List of atomIDs in new_lig_top to add to core_top
            atom_mappings: List [HybridAtoms]
                List of HybridAtoms object keeping track of the atom renumbering
            
        Returns
        --------
            core_top, atom_mappings: Topology, List[HybridAtoms]
                core_top is the updated topology object
                atom_mappings is the list keeping track of the atom renumbering
    """

    # Since this will be called just after adding the bonds, 
    # we assume that the number of residues in the topology equals the residue we are currently working for.
    # NOTE: might have to change this if I make the core different than the -R groups? 
    mres = core_top.SOLUTEATOM.content[-1].MRES
    # Debugging statement to see things more clearly in the file. 
    
    debugging_mode = False
    
    if debugging_mode:
        newAngle = blocks.bondangle_type(IT=999999, JT=999999, KT=999999, ICT=99999)
        core_top.BONDANGLE.content.append(newAngle)
        core_top.BONDANGLE.NTHE += 1
    
    # Here we want to include all bonds that have at least one atom 
    # that is part of the -R group. 
    for angle in new_lig_top.BONDANGLE.content:
        if angle.IT not in new_lig_rgroup and angle.JT not in new_lig_rgroup and \
           angle.KT not in new_lig_rgroup: continue
           
        # note: we need to atom_mappings to be properly set here.
        newAngle = blocks.bondangle_type(IT=find_new_id(atom_mappings, mres, angle.IT),
                                         JT=find_new_id(atom_mappings, mres, angle.JT),
                                         KT=find_new_id(atom_mappings, mres, angle.KT),
                                         ICT=angle.ICT)    
        # Append this new bond to the topology
        core_top.BONDANGLE.append(newAngle)
        core_top.BONDANGLE.NTHE += 1
    
    # Do the same thing for bonds with hydrogen: 
    
    if debugging_mode:
        newAngle = blocks.bondangle_type(IT=999999, JT=999999, KT=999999, ICT=99999)
        core_top.BONDANGLEH.content.append(newAngle)
        core_top.BONDANGLEH.NTHEH += 1
    
    for angle in new_lig_top.BONDANGLEH.content:
        if angle.IT not in new_lig_rgroup and angle.JT not in new_lig_rgroup and \
           angle.KT not in new_lig_rgroup: continue
           
        # note: we need to atom_mappings to be properly set here.
        newAngle = blocks.bondangle_type(IT=find_new_id(atom_mappings, mres, angle.IT),
                                         JT=find_new_id(atom_mappings, mres, angle.JT),
                                         KT=find_new_id(atom_mappings, mres, angle.KT),
                                         ICT=angle.ICT)  
        
        # Append this new bond to the topology
        core_top.BONDANGLEH.append(newAngle)
        core_top.BONDANGLEH.NTHEH += 1
    
    return core_top, atom_mappings

def addHybridTorsions(core_top, new_lig_top, new_lig_rgroup, atom_mappings):
    """
        This function will add all of the angles from a new 
        ligand into the topology with the core. 
        
        This includes the bonds that the -R group forms with itself
        + the one connecting bond between the core region and the R group.
    
        Parameters
        ----------
            core_top: topology
                Core gromos topology object which will be updated (atoms added to this object)
            new_lig_top: topology
                Gromos topology object of the other ligand for which we will add atoms to core_top
            new_lig_rgroup: List [int]
                List of atomIDs in new_lig_top to add to core_top
            atom_mappings: List [HybridAtoms]
                List of HybridAtoms object keeping track of the atom renumbering
            
        Returns
        --------
            core_top, atom_mappings: Topology, List[HybridAtoms]
                core_top is the updated topology object
                atom_mappings is the list keeping track of the atom renumbering
    """

    # Since this will be called just after adding the bonds, 
    # we assume that the number of residues in the topology equals the residue we are currently working for.
    # NOTE: might have to change this if I make the core different than the -R groups? 
    mres = core_top.SOLUTEATOM.content[-1].MRES
    # Debugging statement to see things more clearly in the file. 
    
    debugging_mode = False
    
    if debugging_mode:
        newTorsion = blocks.top_dihedral_type(IP=99999, JP=99999, KP=99999, LP=99999, ICP=99999)
        core_top.DIHEDRAL.content.append(newTorsion)
        core_top.DIHEDRAL.NPHI += 1
    
    # Here we want to include all torsions that have at least one atom 
    # that is part of the -R group. 
    for torsion in new_lig_top.DIHEDRAL.content:
        if torsion.IP not in new_lig_rgroup and torsion.JP not in new_lig_rgroup and \
           torsion.KP not in new_lig_rgroup and torsion.LP not in new_lig_rgroup: continue
           
        # note: we need to atom_mappings to be properly set here.
        newTorsion = blocks.top_dihedral_type(IP=find_new_id(atom_mappings, mres, torsion.IP),
                                              JP=find_new_id(atom_mappings, mres, torsion.JP),
                                              KP=find_new_id(atom_mappings, mres, torsion.KP),
                                              LP=find_new_id(atom_mappings, mres, torsion.LP),
                                              ICP=torsion.ICP)    
        # Append this new bond to the topology
        core_top.DIHEDRAL.append(newTorsion)
        core_top.DIHEDRAL.NPHI += 1
    
    # Do the same thing for torsions with hydrogen: 

    if debugging_mode:
        newTorsion = blocks.dihedralh_type(IPH=99999, JPH=99999, KPH=99999, LPH=99999, ICPH=99999)
        core_top.DIHEDRALH.content.append(newTorsion)
        core_top.DIHEDRALH.NPHIH += 1
    
    for torsion in new_lig_top.DIHEDRALH.content:
        if torsion.IPH not in new_lig_rgroup and torsion.JPH not in new_lig_rgroup and \
           torsion.KPH not in new_lig_rgroup and torsion.LPH not in new_lig_rgroup: continue
           
        # note: we need to atom_mappings to be properly set here.
        newTorsion = blocks.dihedralh_type(IPH=find_new_id(atom_mappings, mres, torsion.IPH),
                                           JPH=find_new_id(atom_mappings, mres, torsion.JPH),
                                           KPH=find_new_id(atom_mappings, mres, torsion.KPH),
                                           LPH=find_new_id(atom_mappings, mres, torsion.LPH),
                                           ICPH=torsion.ICPH)    
        # Append this new bond to the topology
        core_top.DIHEDRALH.append(newTorsion)
        core_top.DIHEDRALH.NPHIH += 1
    
    return core_top, atom_mappings


def addHybridImproperTorsions(core_top, new_lig_top, new_lig_rgroup, atom_mappings):
    """
        This function will add all of the angles from a new 
        ligand into the topology with the core. 
        
        This includes the bonds that the -R group forms with itself
        + the one connecting bond between the core region and the R group.
    
        Parameters
        ----------
            core_top: topology
                Core gromos topology object which will be updated (atoms added to this object)
            new_lig_top: topology
                Gromos topology object of the other ligand for which we will add atoms to core_top
            new_lig_rgroup: List [int]
                List of atomIDs in new_lig_top to add to core_top
            atom_mappings: List [HybridAtoms]
                List of HybridAtoms object keeping track of the atom renumbering
            
        Returns
        --------
            core_top, atom_mappings: Topology, List[HybridAtoms]
                core_top is the updated topology object
                atom_mappings is the list keeping track of the atom renumbering
    """

    # Since this will be called just after adding the bonds, 
    # we assume that the number of residues in the topology equals the residue we are currently working for.
    # NOTE: might have to change this if I make the core different than the -R groups? 
    mres = core_top.SOLUTEATOM.content[-1].MRES
    # Debugging statement to see things more clearly in the file. 
    
    debugging_mode = False
    
    if debugging_mode:
        newTorsion = blocks.impdihedral_type(IQ=99999, JQ=99999, KQ=99999, LQ=99999, ICQ=99999)
        core_top.IMPDIHEDRAL.content.append(newTorsion)
        core_top.IMPDIHEDRAL.NQHI += 1
    
    # Here we want to include all improper torsions that have at least one atom 
    # that is part of the -R group. 
    for torsion in new_lig_top.IMPDIHEDRAL.content:
        if torsion.IQ not in new_lig_rgroup and torsion.JQ not in new_lig_rgroup and \
           torsion.KQ not in new_lig_rgroup and torsion.LQ not in new_lig_rgroup: continue
           
        # note: we need to atom_mappings to be properly set here.
        newTorsion = blocks.impdihedral_type(IQ=find_new_id(atom_mappings, mres, torsion.IQ),
                                             JQ=find_new_id(atom_mappings, mres, torsion.JQ),
                                             KQ=find_new_id(atom_mappings, mres, torsion.KQ),
                                             LQ=find_new_id(atom_mappings, mres, torsion.LQ),
                                             ICQ=torsion.ICQ)    
        # Append this new bond to the topology
        core_top.IMPDIHEDRAL.append(newTorsion)
        core_top.IMPDIHEDRAL.NQHI += 1
    
    # Do the same thing for torsions with hydrogen: 

    if debugging_mode:
        newTorsion = blocks.impdihedralh_type(IQH=99999, JQH=99999, KQH=99999, LQH=99999, ICQH=99999)
        core_top.IMPDIHEDRALH.content.append(newTorsion)
        core_top.IMPDIHEDRALH.NQHIH += 1
    
    for torsion in new_lig_top.IMPDIHEDRALH.content:
        if torsion.IQH not in new_lig_rgroup and torsion.JQH not in new_lig_rgroup and \
           torsion.KQH not in new_lig_rgroup and torsion.LQH not in new_lig_rgroup: continue
           
        # note: we need to atom_mappings to be properly set here.
        newTorsion = blocks.impdihedralh_type(IQH=find_new_id(atom_mappings, mres, torsion.IQH),
                                           JQH=find_new_id(atom_mappings, mres, torsion.JQH),
                                           KQH=find_new_id(atom_mappings, mres, torsion.KQH),
                                           LQH=find_new_id(atom_mappings, mres, torsion.LQH),
                                           ICQH=torsion.ICPH)    
        # Append this new bond to the topology
        core_top.IMPDIHEDRALH.append(newTorsion)
        core_top.IMPDIHEDRALH.NQHIH += 1
    
    return core_top, atom_mappings

def addLigandToTopology(core_top, new_lig_top, new_lig_rgroup, atom_mappings):
    """
    This function will add all of the atoms from a new 
    ligand into the topology with the core. 
    
    Returns a new topology object which contains the hybrid topology
    
    Parameters
    ----------
        core_top: topology contains the first few ligands
        
        new_lig_top: topology of the new ligand for which we will 
                     add some parameters
           
        new_lig_rgroup: List [int]
            indices of -R group atoms of core molecule 
    """
    
    # Add Bonds, Angles, etc.
   
    (core_top, atom_mappings) = addHybridAtoms(core_top, new_lig_top, new_lig_rgroup, atom_mappings)
    (core_top, atom_mappings) = addHybridBonds(core_top, new_lig_top, new_lig_rgroup, atom_mappings)
    (core_top, atom_mappings) = addHybridAngles(core_top, new_lig_top, new_lig_rgroup, atom_mappings)
    (core_top, atom_mappings) = addHybridTorsions(core_top, new_lig_top, new_lig_rgroup, atom_mappings)
    (core_top, atom_mappings) = addHybridImproperTorsions(core_top, new_lig_top, new_lig_rgroup, atom_mappings)
    
    return core_top, atom_mappings

#
# Helper functions
#

def get_bond_list(topology):
    """
        returns a list of tuples 
        corresoponding to all bonds 
        in the topology. 
    """
    bonds = [] # will be a list of tuples
    
    for bond in topology.BOND.content:
        bonds.append( (bond.IB, bond.JB) )
    for bond in topology.BONDH.content:
        bonds.append( (bond.IB, bond.JB) )
    
    return bonds

def find_12_neighbours(atomID, bonds):
    """
        This will find all of the directly bonded atoms to atomID
        
        Returns
        ------- 
            neigh12: List [int]
                list of neighbours
    """
    neigh12 = [] 
    for a, b in bonds:
        if a == atomID: neigh12.append(b)
        elif b == atomID: neigh12.append(a)
    
    return sorted(list(set(neigh12)))


def find_13_neighbours(atomID, bonds, neigh12):
    """
        This will find all of 1,3 neighbours of atomID
        
        Returns
        ------- 
            neigh13: List [int]
                list of 1,3 neighbours
    """

    neigh13 = []      
    for a, b in bonds:
        # skip all bonds which don't include a member of neigh12
        if a not in neigh12 and b not in neigh12: continue
        if a == atomID or b == atomID: continue # don't add the atom itself
        
        # If we reach here, figure out which one is not the second neighbour
        if a in neigh12: neigh13.append(b)
        else: neigh13.append(a)
    
    # Remove potential duplicates, and sort
    return sorted(list(set(neigh13)))

def find_14_neighbours(atomID, bonds, neigh12, neigh13):
    """
        This will find all of 1,4 neighbours of atomID
        
        Returns
        ------- 
            neigh14: List [int]
                list of 1,4 neighbours
    """
    neigh14 = []
    for a, b in bonds:
        # skip all bonds which don't include a member of neigh13
        if a not in neigh13 and b not in neigh13: continue
        if a in neigh12 or b in neigh12: continue
        
        # If we reach here, figure out which one is a third neighbour
        if a in neigh13: neigh14.append(b)
        else: neigh14.append(a)
        
    # Now we need to remove any potential "false" 3rd neighbour
    # which occur for cyclic systems.
    
    for neigh in neigh14:
        if neigh in neigh13 or neigh in neigh12: neigh14.remove(neigh)
    
    return sorted(list(set(neigh14)))


def findLastAtomsOfResidues(top):
    """
    This function finds the last atom of every residue from the topology, 
    and returns this as a list of ints (atom IDs).
    """  
    last_atoms = [] # list to append results to
    prev_res = 1
    for atom in top.SOLUTEATOM.content:
        if atom.MRES != prev_res: 
            last_atoms.append(atom.ATNM-1)
            prev_res = atom.MRES
    
    # Always need to append the last one manually
    last_atoms.append(top.SOLUTEATOM.content[-1].ATNM)
    return last_atoms

def addExclusions(core_top, lig1_atoms):
    """
        This function will add back the exclusions (1,2/1,3 intra/inter ligand)
        as well as the 1,4 between the core of molecule 1 and the new side chains that 
        were added to the molecule.
        
        Before we start writing this function we have:
        
        Ligand1 core with all its exclusions with itself ok.
            --> Add the 1,2/1,3 and 1,4 with all additional ligand rgroups
        
        Ligand1 rgroup with all its exclusions with itself ok.
            --> Add the exclusions with all atoms of the new R groups that were added (INE)
            --> all of its INE14 are already properly setup because they are with itself and with its own core.
            
        Ligands 2 ... N:
            --> Add the full exclusions with other ligands (with larger atomID)
            --> Add the INE14s with itself (not with the core, with itself)
                                            (be careful to exclude fake 1,4 of other hybrid parts)
            
    """
    
    core_atoms, lig_r1_atoms = lig1_atoms    
    bonds = get_bond_list(core_top)
    
    last_atoms = findLastAtomsOfResidues(core_top)
    
    for atom in core_top.SOLUTEATOM.content:
        # find neighbours:
        neigh12 = find_12_neighbours(atom.ATNM, bonds)        
        neigh13 = find_13_neighbours(atom.ATNM, bonds, neigh12)
        neigh14 = find_14_neighbours(atom.ATNM, bonds, neigh12, neigh13)
        
        # case 1 - ligand 1 core
        if atom.MRES == 1 and atom.ATNM in core_atoms:
            
            # Here we will just recalculate everything
            atom.INEvalues = sorted([x for x in neigh12+neigh13 if x > atom.ATNM])
            atom.INE = len(atom.INEvalues)
            
            atom.INE14values = sorted([x for x in neigh14 if x > atom.ATNM])
            atom.INE14 = len(atom.INE14values)
        
        else: # -R groups
            
            # We need to add to the 1,2/1,3 exclusions all atoms belonging to residues N and above 
            full_excl = sorted(neigh12+neigh13)
            
            # Do not add the exclusions with other ligands to test (see if ptp does that automatically)
            # to do so, uncomment next line
            #full_excl.extend(range(last_atoms[atom.MRES-1]+1, last_atoms[-1]+1))
            
            atom.INEvalues = [x for x in full_excl if x > atom.ATNM]
            atom.INE = len(atom.INEvalues)
            
            # For the 1,4 exclusions, include only what is part of the same ligand. 
            # 1,4 exclusions with the core for -R of lig1 == same residue 
            # 1,4 exclusions with the core for -R of lig2...N == included
            # in the INE14 for the core which are present first in the topology.
            
            bounds = (atom.ATNM+1, last_atoms[atom.MRES-1])
            
            neigh14_updated = [x for x in neigh14  if x >= bounds[0] and x <= bounds[1]]
                
            atom.INE14values = sorted(neigh14_updated)
            atom.INE14 = len(atom.INE14values)

            
    return core_top

#
# Smaller functions
#

def adjustMasses(topology):
    """
        This function will adjust the masses of ligands, 
        to ensure they match gromos convention.
    
        This was a problem with openFF ligands, which had slightly different
        masses, which protein part and ligand part did not have the same 
        masses.
    
    """
    
    for atom in topology.SOLUTEATOM.content:
        if   atom.MASS == 1.00795 : atom.MASS = 1.008
        elif atom.MASS == 12.01078: atom.MASS = 12.01
        elif atom.MASS == 14.00672: atom.MASS = 14.01
        elif atom.MASS == 15.99943: atom.MASS = 16.0
        
        # Maybe add chlorines, fluorines, etc. later
            
    return topology

#
# Main functions coordinating everything
#

def constructPerturbedTopology(hybrid_topology, lig1_core, out_path, alphalj = 1.0, alphacrf = 1.0):
    """
        This function will create the perturbed topology for our hybrid topology. 
        All ligands atoms will be written down (so energies in the EDS blocks of the output match what we want).
        
        Shared core atoms will have the same atom type code in all states. 
        -R groups will be dummy in all other states. 
    
        note: The core isn't sorted so there are ligand 1 -R group atoms in the middle of the ptp.
        
        Parameters
        ----------
            hybrid_topology: pyGromos Top 
                regular topology from which we will create the hybrid topology
            lig1_core: List [int]
                list of atoms which make up the core (present for all perturbed states)
            out_path: str
                path to save the output ptp in
                
            alphalj: float
                alpha parameter for the LJ interaction
            alphacrf: float
                    alpha parameter for the CRF interaction
    
    """    
    # count number of residues
    num_atoms = hybrid_topology.SOLUTEATOM.content[-1].ATNM
    num_states = hybrid_topology.SOLUTEATOM.content[-1].MRES
    
    dummy_iac = hybrid_topology.ATOMTYPENAME.content[0][0] # assumes dummy iac is always last
    
    # open the output file.
    f = open(out_path, 'w')

    # Write title block.
    
    date_time = datetime.now().strftime("%d/%m/%Y %H:%M")    
    f.write('TITLE\n\tFile created automatically from ' + hybrid_topology.path +'\n')
    f.write('\tby ' + os.environ['USER'] + ' on ' + date_time + '\nEND\n')
    
    # Write the MPERTATOM block
    
    f.write('MPERTATOM\n# NJLA:   number of perturbed atoms\n')
    f.write('# NPTB:   number of listed perturbation (i.e. number of perturbation states)\n')
    f.write('# NJLA    NPTB\n')
    
    # Write the values to NJLA and NPTB:
    f.write('\t' + str(num_atoms) + '\t\t' + str(num_states) + '\n')
    
    # Write state identifiers
    identifiers = ""
    for i in range(num_states) : identifiers += "ligand"+str(i+1) + "\t"
    
    f.write('# identifiers of the states\n\t' + identifiers +'\n')
    
    # comment to understand the content
    f.write('#  NR  NAME IAC(1) CHARGE(1) ...  IAC(n) CHARGE(n)  ALPHLJ  ALPHCRF\n')
    
    # Loop over topology atoms to include in the ptp
    
    for atom in hybrid_topology.SOLUTEATOM.content:
        
        ptp_str = '\t' + str(atom.ATNM) + '\t' + str(atom.PANM) +'\t'
        
        # charge offset to make formating look good
        q_offset = ' ' if atom.CG > 0 else ''
        
        if atom.MRES == 1 and atom.ATNM in lig1_core: 
            for i in range(num_states): ptp_str += str(atom.IAC) + '\t' + q_offset + "{:.5f}".format(atom.CG) + '\t'

        else: # Dealing with -R groups
            for i in range(num_states):
                if i == atom.MRES-1:
                    ptp_str += str(atom.IAC) + '\t' + q_offset + "{:.5f}".format(atom.CG) + '\t'
                else:
                    ptp_str += str(dummy_iac) + '\t ' + "{:.5f}".format(0) + '\t'
        
        # Always prepend alpha LJ and alpha CRF values
        
        ptp_str += "{:.5f}".format(alphalj) + '\t' + "{:.5f}".format(alphacrf) + '\n'
        f.write(ptp_str)
        
    # close the block
    f.write('END\n')
    f.close()
    
    return None
    
def constructHybridConformation(new_ligand_tops, connecting_points, paths_input_cnfs, path_out_cnf):
    """
        This function will create a Hybrid Conformation, based on the hybrid topology, 
        and single ligand conformations given. 
        
        Note: This functions needs to be adapted so we can make a hybrid conformation for the 
              complex!
        
        Parameters
        ----------
            new_ligand_tops: List [str]
                path to the individual topologies (of the R groups)
            connecting_points: List [ Tuple(List[int],List[int])]    
                List of the "connecting point" between core and -R group
                The first inner lists correspond to the atom IDs of the atoms belonging 
                to the core directly bound to a -R group atom.
                The second inner list corresponds to the atom IDs of the atoms belonging
                to the -R group directly bound to the core. 
            paths_input_cnfs: List [str]
                list of paths of the input conformations to combine (pre-aligned ideally)
            path_out_cnf: str
                path of the output cnf
        Returns
        --------
            None
        
    """

    core_cnf = Cnf(paths_input_cnfs[0])
    
    new_hybrid_cnf = copy.deepcopy(core_cnf)
    
    # Remove all waters (We will add them back at the end) 
    new_hybrid_cnf.delete_residue(resName="WAT")
    new_hybrid_cnf.delete_residue(resName="SOLV")
    
    new_hybrid_pos = new_hybrid_cnf.POSITION.content
    
    curAtomID = new_hybrid_cnf.POSITION.content[-1].atomID +1
    curResID  = new_hybrid_cnf.POSITION.content[-1].resID +1
    
    # Add all additional -R groups
    for i, cnf_path in enumerate(paths_input_cnfs[1:]):
        newResId = i+2
        lig_cnf =  Cnf(cnf_path)
        
        lig_core, lig_rgroup = findAtomsInCoreAndRGroup(new_ligand_tops[i], connecting_points[i+1])
        
        for pos in lig_cnf.POSITION.content:
            if pos.resID == 2: break
            
            # Otherwise check if atomID matches the -R groups
            if pos.atomID in lig_rgroup:
                # Append a new atom position
                tmp = atomP(resID = curResID, resName = pos.resName, atomType=pos.atomType, 
                            atomID = curAtomID, xp = pos.xp, yp = pos.yp, zp= pos.zp)
        
                new_hybrid_pos.append(tmp)
                curAtomID += 1
        # end of the addition of the ligand
        curResID += 1
        
    # Now that all ligands have been added, add all the water molecules back.
    # They will be renamed SOLV to match gromos conventions.
        
    for pos in core_cnf.POSITION.content:
        if pos.resID == 1: continue 
        
        tmp = atomP(resID = curResID, resName = "SOLV", atomType=pos.atomType, 
                    atomID = curAtomID, xp = pos.xp, yp = pos.yp, zp= pos.zp)
        
        new_hybrid_pos.append(tmp)
        curAtomID += 1
        if pos.atomType == "H2": curResID += 1 # Each time we have a new water
    
    new_hybrid_cnf.write(path_out_cnf)
    new_hybrid_cnf.write_pdb(path_out_cnf.replace('.cnf', '.pdb'))
    return None


def constructHybridTopology(core_top, new_ligand_tops, atom_mappings, connecting_points, path_out_top):
    """
    This is the main function called to construct Hybrid topologies, 
    which will call all subfunctions doing the job.
    
    The user needs to provide the topologies, as well as a manually curated
    list of atom mappings (so the program will figure out how to add the torsions)
    between core and new -R groups.
    
    This mapping procedure may be automated in the future.
    
    Parameters
    ----------
        core_top: pygromos Top
            Topology of the core
        new_ligand_tops: List [Top]
            Topologies of the new ligands which will have the added -R groups
        atom_mappings: List [PerturbedAtoms]
            List keeping track of the mapping of perturbed atoms between their original topology 
            and the new hybrid one.
            
        connecting_points: List [ Tuple(List[int],List[int])]    
            List of the "connecting point" between core and -R group
            The first inner lists correspond to the atom IDs of the atoms belonging 
            to the core directly bound to a -R group atom.
            The second inner list corresponds to the atom IDs of the atoms belonging
            to the -R group directly bound to the core. 
            
        path_out_top: str
            path to save the final topology to.
           
    Returns
    -------
        None
    
    """
    
    new_core = copy.deepcopy(core_top)
    lig1_atoms = findAtomsInCoreAndRGroup(core_top, connecting_points[0])
    
    # 1: Add all of the topologies together
    
    for i, lig_top in enumerate(new_ligand_tops):
        lig_core, lig_rgroup = findAtomsInCoreAndRGroup(lig_top, connecting_points[i+1])     
        print ('Working on addition of ligand ' + str(i+2) + ' has: ' + str(len(lig_rgroup)) + ' atoms.')
        (new_core, atom_mappings) = addLigandToTopology(new_core, lig_top, lig_rgroup, atom_mappings)        
        
    # 2: Reset the proper exclusions given new bonds added.
    new_core = addExclusions(new_core, lig1_atoms)
    
    # 3: Adjust small things in some blocks
    
    # Change masses (openFF assigns very slightly different masses)
    new_core = adjustMasses(new_core)

    # Last atom needs to be the end of a charge group.
    new_core.SOLUTEATOM.content[-1].CGC = 1
    
    num_atoms = new_core.SOLUTEATOM.content[-1].ATNM
    
    # solute-molecules, pressure and temperature groups 
    # need to be readjusted. 
    new_core.SOLUTEMOLECULES.content = [['1'], [str(num_atoms)]]
    
    new_core.PRESSUREGROUPS.content = [['1'], [str(num_atoms)]]
    new_core.TEMPERATUREGROUPS.content = [['1'], [str(num_atoms)]]
    
    # Adjust residue names:
    
    new_core.RESNAME.content = []
    new_core.RESNAME.content.append([str(new_core.SOLUTEATOM.content[-1].MRES)])
    for i in range(1, new_core.SOLUTEATOM.content[-1].MRES+1): 
        new_core.RESNAME.content.append(['LI'+str(i)])
    
    # Adjust TITLE block
    
    date_time = datetime.now().strftime("%d/%m/%Y %H:%M")    
    new_core.TITLE.content = ['\tFile created automatically by ' + os.environ['USER'] + ' on ' + date_time]
    
    # Save results
    
    new_core.write(path_out_top)
    
    # 4: Make a Perturbed topology which matches this.
    
    perturbed_top_path = new_core.path.replace(".top", ".ptp")
    constructPerturbedTopology(new_core, lig1_atoms[0], perturbed_top_path)
    
    return None













