#!/usr/bin/env python
"""
Utility for selecting molecules from a .xyz file

The default output is called out.xyz but can be specified with -o

Usage:
pick_mol.py input_file.xyz 7 13 42

"""
import sys
import argparse
from copy import copy

from fromage.io import read_file as rf
from fromage.io import edit_file as ef


def picker(in_name, out_name, labels, bonding, thresh, reverse=False):
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
    atoms = rf.mol_from_file(in_name)
    if thresh == 999:
        thresh = None
    atoms.set_bonding(bonding=bonding, thresh=thresh)

    # if the label does not exist
    if max(labels) > len(atoms):
        raise ValueError("One or more atom labels were too high!")

    selected = atoms.copy()
    selected.atoms = []

    for label in labels:
        mol = atoms.select(label)

        # prevent double selecting
        if mol[0] in selected:
            raise ValueError("Atom " + str(label) + " was already selected!")

        selected += mol

    # to print out the non specified atoms
    if reverse:
        to_remove = [copy(i) for i in selected]
        selected = atoms.copy()
        selected.atoms = []
        for atom in atoms:
            if atom not in to_remove:
                selected.append(atom)

    selected.write_xyz(out_name)

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input .xyz file",
                        default="cluster.xyz")
    parser.add_argument("-o", "--output", help="Name of the output file",
                        default="out.xyz", type=str)
    parser.add_argument("labels", help="Number label(s) for the atom(s) in the picked molecule(s)",
                        default=[0], type=int, nargs='*')
    parser.add_argument(
        "-b", "--bonding", help="Type of bonding used for detecting full molecules. The options are dis, cov and vdw", default="dis", type=str)
    parser.add_argument("-T", "--thresh", help="Threshold distance to select a bond. The default pairs are dis:1.8, cov:0.1, vdw:0.3. To select default, just use a threshold of 999", default=999, type=float)
    parser.add_argument(
        "-r", "--reverse", help="Print all atoms except selected molecules", action="store_true")
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)

    # users will use atom labels starting from 1
    new_labels = [a - 1 for a in args.labels]

    # call the main function
    picker(args.input, args.output, new_labels,
           args.bonding, args.thresh, reverse=args.reverse)
