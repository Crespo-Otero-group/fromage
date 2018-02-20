#!/usr/bin/python3
"""
Utility for selecting molecules from a .xyz file

The default output is called out.xyz but can be specified with -o

Usage:
pick_mol.py input_file.xyz 7 13 42
"""
import sys
import argparse

from cryspy.io import read_file as rf
from cryspy.io import edit_file as ef
from cryspy.utils import handle_atoms as ha


def picker(in_name, out_name, labels, max_bl):
    """
    Pick out molecules from an xyz file

    Parameters
    ----------
    in_name : str
        Name of the .xyz input file
    out_name : str
        Name of the .xyz output file
    labels : list of ints
        Labels of the atoms in the molecules to select (array starts at 0)
    max_bl : float
        Maximun distance in Angstrom that constitutes a bond
    """
    atoms = rf.read_xyz(in_name)[-1]

    # if the label does not exist
    if max(labels) > len(atoms):
        raise ValueError("One or more atom labels were too high!")

    selected = []
    for label in labels:
        mol = ha.select(max_bl, atoms, label)

        # prevent double selecting
        if mol[0] in selected:
            raise ValueError("Atom " + str(label) + " was already selected!")

        selected += mol
    ef.write_xyz(out_name, selected)

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input .xyz file",
                        default="cluster.xyz")
    parser.add_argument("-o", "--output", help="Name of the output file",
                        default="out.xyz", type=str)
    parser.add_argument("labels", help="Number label(s) for the atom(s) in the picked molecule(s)",
                        default=[0], type=int, nargs='*')
    parser.add_argument("-b", "--bond", help="Maximum length in Angstrom that qualifies as a bond",
                        default=1.6, type=float)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)

    # users will use atom labels starting from 1
    new_labels = [a - 1 for a in args.labels]

    # call the main function
    picker(args.input, args.output, new_labels, args.bond)
