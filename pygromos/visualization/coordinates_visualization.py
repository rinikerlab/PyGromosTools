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

    if len(protein_resi) > 0 or len(water_resi) > 0 or len(ions_resn) > 0:
        representation = [
            {"type": "cartoon", "params": {"sele": " ".join(map(str, protein_resi)), "color": "residueindex"}},
        ]
        view = nj.show_mdtraj(traj, representations=representation)

    else:
        view = nj.show_mdtraj(traj)

    if len(protein_resn) > 0:
        view.add_representation("cartoon", selection=" ".join(map(str, protein_resn)), color="residueindex")

        if len(protein_resn) < 16:
            view.add_representation("hyperball", selection=" ".join(map(str, protein_resn)))
    if len(water_resi) > 0 and traj.n_frames < 300:
        view.add_representation("line", selection=" ".join(map(str, water_resi)))
    if len(ions_resn) > 0:
        view.add_representation("hyperball", selection=" ".join(map(str, ions_resn)), radius=3, color="green")
    return view
