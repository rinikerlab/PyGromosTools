def generate_dual_topology_approach(cnfA, cnfB, topA, topB, eds:bool=False):
    ##Atom Mapping
    atom_mappingAB, smart = find_atom_mapping(cnfA=cnfA, cnfB=cnfB)
    
    ##Coordinates
    cnfA, cnfB = align_cnfs_with_MCS(cnfA=cnfA, cnfB=cnfB, atom_mappingAB=atom_mappingAB)
    cnf_comb = copy.deepcopy(cnfA)
    #cnf_comb += cnfB # needs to be implemented
    
    for pos in cnfB.POSITION:
        cnf_comb.POSITION.append(pos)
    
    cnf_comb.supress_atomPosition_singulrarities()


    ##Top
    nAtoms_top1 = len(top1.SOLUTEATOM)
    top_comb = copy.deepcopy(top1)
    top_comb += top2
    
    ### Pertubation
    ptp_comb = ptp.Pertubation_topology()
    
    if(eds):
        from pygromos.files.blocks.pertubation_blocks import MPERTATOM
        from pygromos.files.blocks.pertubation_blocks import  atom_eds_pertubation_state, pertubation_eds_state

        tops = [topA, topB]
        dummyState = pertubation_eds_state(IAC=22, CHARGE=0)

        numStates=len(tops)
        IND = 1
        atom_states = []
        for top_ind, top in enumerate(tops):
            for atom in top.SOLUTEATOM:        
                states = {}
                for ctop in range(ntops):
                    if(ctop==top_ind):
                        states.update({ctop+1:pertubation_eds_state(IAC=atom.IAC, CHARGE=atom.CG)})
                    else:
                        states.update({ctop+1:dummyState})

                atom_ptp = atom_eds_pertubation_state(NR=IND, NAME=atom.PANM, STATES=states)
                atom_states.append(atom_ptp)
                IND+=1

        ptp_comb.add_block(block=MPERTATOM(NJLA=len(atom_states), NPTB=numStates, STATEATOMS=atom_states))

    else:
        from pygromos.files.blocks.pertubation_blocks import PERTATOMPARAM
        from pygromos.files.blocks.pertubation_blocks import  atom_lam_pertubation_state, pertubation_lam_state_nonbonded

        tops = [topA, topB]
        build_dummyState = lambda m: pertubation_lam_state_nonbonded(IAC=22, CHARGE=0, MASS=m)

        numStates=len(tops)
        IND = 1
        atom_states = []
        for top_ind, top in enumerate(tops):
            for atom in top.SOLUTEATOM:        
                states = {}
                for ctop in range(ntops):
                    if(ctop==top_ind):
                        states.update({ctop+1:pertubation_lam_state_nonbonded(IAC=atom.IAC, CHARGE=atom.CG, MASS=atom.MASS)})
                    else:
                        states.update({ctop+1:build_dummyState(atom.MASS)})

                atom_ptp = atom_lam_pertubation_state(NR=IND, RES=atom.MRES, NAME=atom.PANM, STATES=states,)
                atom_states.append(atom_ptp)
                IND+=1

        ptp_comb.add_block(block=PERTATOMPARAM(NJLA=len(atom_states), STATEATOMS=atom_states))
    
    return cnf_comb, top_comb, ptp_comb

def generate_hybrid_topology_approach(cnfA, cnfB, topA, topB):
    ##Atom Mapping
    atom_mappingAB, smart = find_atom_mapping(cnfA=cnfA, cnfB=cnfB)
    
    ##Coordinates
    cnfA, cnfB = align_cnfs_with_MCS(cnfA=cnfA, cnfB=cnfB, atom_mappingAB=atom_mappingAB)
    cnf_comb, present_atoms = merge_states(cnfA=cmol1, cnfB=cmol2, atomMatchingAB=atom_mappingAB, dist_tresh=0.0, _doNotChangeAtomType=True, _doUpdateAtomMapping=True) #no distance collapsing
    
    ##Top
    
    
    ### Pertubation
    
    
    return cnf_comb
    #return cnf_comb, top_comb, ptp_comb

def generate_single_topology_approach(cnfA, cnfB, topA, topB):
    ##Atom Mapping
    atom_mappingAB, smart = find_atom_mapping(cnfA=cnfA, cnfB=cnfB)
    
    ##Coordinates
    cnfA, cnfB = align_cnfs_with_MCS(cnfA=cnfA, cnfB=cnfB, atom_mappingAB=atom_mappingAB)
    cnf_comb,present_atoms = merge_states(cnfA=cmol1, cnfB=cmol2, atomMatchingAB=atom_mappingAB, dist_tresh=0.09) #no distance collapsing
    

    return cnf_comb #, top_comb, ptp_comb