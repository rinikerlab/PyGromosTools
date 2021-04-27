"""
    in this file coordinate visualizations are implemented
"""
from pygromos.files.coord.cnf import Cnf
from pygromos.files.trajectory.trc import Trc

import py3Dmol

def show_cnf(cnf: Cnf):
    xyz_str = cnf.get_xyz()
    view=py3Dmol.view(width=400, height=400)
    view.addModel(xyz_str)
    view.setStyle({"stick":{}})
    view.zoomTo()

    return view



def show_coordinate_traj(trc:Trc, cnf: Cnf):
    traj = trc.get_pdb(cnf)
    view = py3Dmol.view(width=400, height=400)
    view.addModelsAsFrames(traj)
    view.setStyle({'model': -1}, {"stick": {}})
    view.animate({"loop": "forwardAndBackward"})

    view.zoomTo()
    return view
