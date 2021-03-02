import unittest
import numpy as np
from itertools import combinations
from pygromos.files.blocks.topology_blocks import atom_pertubation_state, pertubation_state
from pygromos.files.topology.top import Pertubation_topology as ptp

import os
path= os.path.dirname(__file__) +"/testfiles/ptp/eds.ptp"
outpath=os.path.dirname(__file__) +"/testfiles/ptp/out_ptp.ptp"
outpath_less_state=os.path.dirname(__file__) +"/testfiles/ptp/out_ptp_lessStates.ptp"
outpath_less_atom=os.path.dirname(__file__) +"/testfiles/ptp/out_ptp_lessAtoms.ptp"
outpath_new_atom=os.path.dirname(__file__) +"/testfiles/ptp/out_ptp_newAtoms.ptp"
outpath_new_build=os.path.dirname(__file__) +"/testfiles/ptp/out_ptp_newBuild.ptp"
outpath_new_build_complex=os.path.dirname(__file__) +"/testfiles/ptp/out_ptp_newBuild_complex.ptp"

class test_ptp(unittest.TestCase):
    file_class = ptp

    def test_IO(self):
        ptp = self.file_class(path)
        print(ptp)
        ptp.write(outpath)

    def test_get_states(self):
        ptp = self.file_class(path)
        states = ptp.MPERTATOM.states
        print(states)


    def test_add_state(self):
        from pygromos.files.blocks.topology_blocks import  atom_pertubation_state, pertubation_state
        ptp = self.file_class(path)

        on_state = pertubation_state(IAC=16, CHARGE=-1.0)
        new_atoms_state = [atom_pertubation_state(NR=x, NAME="H", STATES={7:on_state, 13:on_state}) if x==1 else atom_pertubation_state(NR=x, NAME="H", STATES={7:on_state})  for x in range(1, 4)]

        ptp.MPERTATOM.add_state_atoms(state_atoms=new_atoms_state)

        ptp.write(outpath_new_atom)


    def test_delete_state(self):
        ptp = self.file_class(path)
        ptp.MPERTATOM.delete_state(stateIDs=6)
        ptp.MPERTATOM.delete_state(stateNames="state3")
        ptp.write(outpath_less_state)

    def test_remove_atoms(self):
        ptp = self.file_class(path)
        ptp.MPERTATOM.delete_atom(atomNR=[3, 4, 5, 6])
        ptp.write(outpath_less_atom)

    def test_new_ptp_from_scratch(self):
        from pygromos.files.blocks.topology_blocks import atom_pertubation_state, pertubation_state
        ptp = self.file_class()
        on_state = pertubation_state(IAC=16, CHARGE=-1.0)

        new_atoms_state = [atom_pertubation_state(NR=x, NAME="H", STATES={x: on_state}) for x in range(1, 4)]

        ptp.MPERTATOM.add_state_atoms(state_atoms=new_atoms_state)
        print(ptp)
        ptp.write(outpath_new_build)

    def test_gen_all_states_ptp(self):
        # INPUT:
        ## the states to use:
        o_state = pertubation_state(IAC=16, CHARGE=-1.0)
        h_state = pertubation_state(IAC=33, CHARGE=-2.0)

        ## Map for molecule ID with Atom IDS
        ### First atom assumed to be O and last two Hs
        molecules_atoms = {1: [1, 2, 3],
                           2: [4, 5, 6],
                           3: [7, 8, 9],
                           4: [10, 11, 12]}

        # BUILD UP STATES
        ## Generate active state mapping:
        max_active_mols_same_time = len(molecules_atoms)
        molecule_states = {}
        state_ID = 1
        for active_mols in range(1, max_active_mols_same_time + 1):
            combis = list(combinations(molecules_atoms, active_mols))
            molecule_states.update({state_ID + ind: x for ind, x in enumerate(combis)})
            state_ID = state_ID + len(combis)
        # gives the state number as key and all the active molecules in this state

        # build state atoms for ptp
        new_atoms_state_dict = {}
        for state, molecules in molecule_states.items():
            for molecule in molecules:
                if (molecule in new_atoms_state_dict):
                    atoms = new_atoms_state_dict[molecule]
                    atoms[0].STATES.update({state: o_state})
                    atoms[1].STATES.update({state: h_state})
                    atoms[2].STATES.update({state: h_state})
                else:
                    atoms = [atom_pertubation_state(NR=molecules_atoms[molecule][0], NAME="O", STATES={state: o_state}),
                             atom_pertubation_state(NR=molecules_atoms[molecule][1], NAME="H1",
                                                    STATES={state: h_state}),
                             atom_pertubation_state(NR=molecules_atoms[molecule][2], NAME="H2",
                                                    STATES={state: h_state})]
                    new_atoms_state_dict.update({molecule: atoms})


        ##finally make a list for our ptp file (#ThanksGromos)
        new_atoms_state = np.concatenate(list(new_atoms_state_dict.values()))

        # BUILD PTP
        ptp = self.file_class(path)
        ptp.MPERTATOM.add_state_atoms(state_atoms=new_atoms_state)
        print(ptp)
        ptp.write(outpath_new_build_complex)

        # TADAAAAA - DONE