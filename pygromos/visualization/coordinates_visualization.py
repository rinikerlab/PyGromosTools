"""
        In this file certain Visualization pretypes could be implemented
"""
import mdtraj
import nglview as nj

from pygromos.utils.amino_acids import three_letter_aa_lib, solvents, ions


def visualize_system(traj: mdtraj.Trajectory) -> nj.NGLWidget:

    protein_resi = set([x.index + 1 for x in traj.top.residues if (x.name in three_letter_aa_lib)])
    protein_resn = set([x.name[:3] for x in traj.top.residues if (x.name in three_letter_aa_lib)])

    water_resi = set([x.index + 1 for x in traj.top.residues if (x.name in solvents)])
    ions_resn = set([x.name for x in traj.top.residues if (x.name in ions)])

    representation = [
        {"type": "cartoon", "params": {"sele": " ".join(map(str, protein_resi)), "color": "residueindex"}},
        {"type": "line", "params": {"sele": " ".join(map(str, water_resi))}},
    ]
    view = nj.show_mdtraj(traj, representations=representation)

    if len(protein_resn) < 16:
        view.add_representation("hyperball", selection=" ".join(map(str, protein_resn)))
    if len(ions_resn) > 0:
        view.add_representation("hyperball", selection=" ".join(map(str, ions_resn)), radius=3, color="green")
    return view
