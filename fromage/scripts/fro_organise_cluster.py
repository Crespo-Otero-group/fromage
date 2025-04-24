#!/usr/bin/env python
"""
Organise a molecular cluster to have the atoms of each molecule
in consecutive lines in the shell.xyz file

Usage:
fro_organise_cluster.py cluster.xyz

Output:
cluster_ordered.xyz

"""
import argparse
import numpy as np
import fromage.io.read_file as rf
from fromage.utils.mol import Mol


def main(cluster_file):
    cluster = rf.mol_from_file(cluster_file)
    list_mols = cluster.segregate(diff_mols = True)

    if not list_mols:
        raise ValueError("No molecules found in the cluster. Check your input files.")

    reordered_atoms = []
    for mol in list_mols:
        for atom in mol:
            reordered_atoms.append(atom)
    new_cluster = Mol(reordered_atoms)
    new_cluster.write_xyz('cluster_ordered.xyz')

if __name__=='__main__':
    parser = ArgumentParser(description='Organise a molecular cluster with molecules in consecutive lines.')
    parser.add_argument('cluster_file', help='Path to the XYZ file of the full cluster')
    args = parser.parse_args()
    main(args.cluster_file)

