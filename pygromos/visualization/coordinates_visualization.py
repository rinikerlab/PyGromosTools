"""
    in this file coordinate visualizations are implemented
"""
from pygromos.files.coord.cnf import Cnf
from pygromos.files.trajectory.trc import Trc
from pygromos.utils.amino_acids import three_letter_aa_lib, ions
from copy import deepcopy
import py3Dmol

def show_cnf(cnf: Cnf):
    cCNF = deepcopy(cnf)

    view=py3Dmol.view(width=400, height=400)

    solute = [resn[:3] for resn in cnf.residues if(cCNF !="SOLV")]
    if(len([res for res in solute if(res in three_letter_aa_lib)])>15 or "SOLV" in cnf.residues):
        solv_cnf = deepcopy(cnf)

        pos = []
        for atomP in cCNF.POSITION:
            atomP.atomType = atomP.atomType[:1]
            if(atomP.resName == "SOLV"):
                pos.append(atomP)
        solv_cnf = deepcopy(cnf)
        solv_cnf.POSITION = pos

        xyz_str = cCNF.get_pdb()
        xyzd_str = cCNF.get_xyz()

        view.addModel(xyz_str)
        view.addModel(xyzd_str)

        view.setStyle({'resn': solute}, {"stick": {}})

        if(len([res for res in solute if(res in three_letter_aa_lib)])>15):
            protein = [res for res in solute if(res in three_letter_aa_lib)]
            view.setStyle({'resn': protein}, {"cartoon": {}})
            view.setStyle({'cartoon': {'arrows': True, 'tubes': True, 'style': 'oval'}})

        view.setStyle({'resn': "SOLV"}, {"line": {}})        # Solvent
        view.setStyle({'resn': ions},  {"sphere": {"color": "lightgreen", "radius":0.7}})        # ions

    else:
        xyz_str = cCNF.get_xyz()
        view.addModel(xyz_str)
        view.setStyle({"model": -1}, {"stick": {}})

    view.zoomTo()
    return view



def show_coordinate_traj(trc:Trc, cnf: Cnf):
    """
    This function visualizes the provided TRC and maps it on the

    Parameters
    ----------
    trc : Trc

    cnf: Cnf

    Returns
    -------

    """
    traj = trc.get_pdb(cnf)
    view = py3Dmol.view(width=400, height=400)
    view.addModelsAsFrames(traj)
    view.setStyle({'model': -1}, {"stick": {}})
    view.animate({"loop": "forwardAndBackward"})

    view.zoomTo()
    return view
