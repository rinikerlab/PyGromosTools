""" Amino Acid Library

Notes
----------
    This module contains all natural amino acids as named tuples.
    Also contained are dictionaries for three letter code and so on.
    Todo: add number of atoms in full atomistic case and other useful properties.

"""
from collections import namedtuple

ions = ["NA+", "CL-"]
solvents = ["SOL", "SOLV", "HOH", "H2O"]
# add here new amino acid features.
amino_acid: namedtuple = namedtuple(
    "amino_acid", ["name", "oneLetter", "threeLetter", "numUnitedAtoms", "numFullAtomistic"]
)

# Add here new Amino acids:
Alanine = amino_acid(name="alanine", oneLetter="A", threeLetter="ALA", numUnitedAtoms=6, numFullAtomistic=None)

Cystein = amino_acid(name="cystein", oneLetter="C", threeLetter="CYSH", numUnitedAtoms=8, numFullAtomistic=None)

AsparticAcid = amino_acid(
    name="aspartic acid", oneLetter="D", threeLetter="ASP", numUnitedAtoms=9, numFullAtomistic=None
)

GlutamicAcid = amino_acid(
    name="glutamic acid", oneLetter="E", threeLetter="GLU", numUnitedAtoms=10, numFullAtomistic=None
)

Phenylalanine = amino_acid(
    name="phenylalanine", oneLetter="F", threeLetter="PHE", numUnitedAtoms=17, numFullAtomistic=None
)

Glycine = amino_acid(name="glycine", oneLetter="G", threeLetter="GLY", numUnitedAtoms=5, numFullAtomistic=None)

Histidine_Neutral = amino_acid(
    name="histidine", oneLetter="H", threeLetter="HIS", numUnitedAtoms=10, numFullAtomistic=None
)
Histidine_HID = amino_acid(
    name="histidine_HID", oneLetter="H", threeLetter="HISA", numUnitedAtoms=13, numFullAtomistic=None
)
Histidine_HIE = amino_acid(
    name="histidine_HIE", oneLetter="H", threeLetter="HISB", numUnitedAtoms=14, numFullAtomistic=None
)
Histidine_positive = amino_acid(
    name="histidine_positive", oneLetter="H", threeLetter="HISP", numUnitedAtoms=15, numFullAtomistic=None
)

Isoleucine = amino_acid(name="isoleucine", oneLetter="I", threeLetter="ILE", numUnitedAtoms=9, numFullAtomistic=None)

Lysine = amino_acid(name="lysine", oneLetter="K", threeLetter="LYS", numUnitedAtoms=12, numFullAtomistic=None)
Lysine_positive = amino_acid(
    name="lysine_positive", oneLetter="K", threeLetter="LYSH", numUnitedAtoms=13, numFullAtomistic=None
)

Leucine = amino_acid(name="leucine", oneLetter="L", threeLetter="LEU", numUnitedAtoms=9, numFullAtomistic=None)

Methionine = amino_acid(name="methionine", oneLetter="M", threeLetter="MET", numUnitedAtoms=9, numFullAtomistic=None)

Asparagine = amino_acid(name="asparagine", oneLetter="N", threeLetter="ASN", numUnitedAtoms=11, numFullAtomistic=None)

Proline = amino_acid(name="proline", oneLetter="P", threeLetter="PRO", numUnitedAtoms=7, numFullAtomistic=None)

Glutamine = amino_acid(name="glutamine", oneLetter="Q", threeLetter="GLN", numUnitedAtoms=12, numFullAtomistic=None)

Arginine = amino_acid(name="arginine", oneLetter="R", threeLetter="ARG", numUnitedAtoms=17, numFullAtomistic=None)

Serine = amino_acid(name="serine", oneLetter="S", threeLetter="SER", numUnitedAtoms=8, numFullAtomistic=None)

Threonine = amino_acid(name="threonine", oneLetter="T", threeLetter="THR", numUnitedAtoms=9, numFullAtomistic=None)

Valine = amino_acid(name="valine", oneLetter="V", threeLetter="VAL", numUnitedAtoms=8, numFullAtomistic=None)

Tryptophane = amino_acid(name="tryptophane", oneLetter="Y", threeLetter="TRP", numUnitedAtoms=21, numFullAtomistic=None)

Tyrosine = amino_acid(name="tyrosine", oneLetter="Y", threeLetter="TYR", numUnitedAtoms=18, numFullAtomistic=None)


# Add new amino Acid into set. rest automatic
aa_set = {
    Alanine,
    Cystein,
    AsparticAcid,
    GlutamicAcid,
    Phenylalanine,
    Glycine,
    Histidine_Neutral,
    Isoleucine,
    Lysine,
    Leucine,
    Methionine,
    Asparagine,
    Proline,
    Glutamine,
    Arginine,
    Serine,
    Threonine,
    Valine,
    Tryptophane,
    Tyrosine,
    Lysine_positive,
    Histidine_HIE,
    Histidine_HID,
    Histidine_positive,
}

# nici
# automatic lib gen
three_letter_aa_lib = {aa.threeLetter: aa for aa in aa_set}
one_letter_aa_lib = {aa.oneLetter: aa for aa in aa_set}
aa_lib = {aa.name: aa for aa in aa_set}
