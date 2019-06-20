#!/usr/bin/env python

"""Functions useful for manipulating atom connectivity.

The original use of these functions is to allow for assign_charges to work and
be able to give a molecule, cluster or periodic cell of atoms the same charges as another
molecule, cluster or periodic cell.

Includes a utility to read a molecule xyz file, a population analysis in g09 and
a target cluster xyz to assign the charges to the cluster.

Usage:
assign_charges.py mol.log clust.xyz

Includes options for Mulliken or RESP and ouptut file names.
"""

import numpy as np
import sys
import argparse

import fromage.io.read_file as rf


def detect_1_connect(in_atoms):
    """
    Make a matrix of first connectivities of a list of atoms.

    Parameters
    ----------
    in_atoms : Mol object
        Atoms which need their connectivity detected
    vectors : 3x3 array-like or None
        Lattice vectors of the system if it is periodic. If not use None

    Returns
    -------
    cnct : numpy matrix
        The matrix where each row and each column correspond to one atom. If the
        two atoms are bonded or the same, the matrix element is 1. Otherwise it
        is 0

    """
    nat_mol = len(in_atoms)
    cnct = np.zeros((nat_mol, nat_mol),dtype=int)
    for i, i_atom in enumerate(in_atoms):
        for j, j_atom in enumerate(in_atoms):
            if np.count_nonzero(in_atoms.vectors) == 0:
                if in_atoms.bonded(i_atom, j_atom):
                    cnct[i][j] = 1
            else:
                if in_atoms.per_bonded(i_atom, j_atom):
                    cnct[i][j] = 1
    return cnct


def expand_connect(in_mat):
    """
    Expand a connectivity matrix

    For one atom, checks which other atoms are connected to it (connectors) and
    which are not (dangles). Then for each dangle checks if it has connectors in
    common with the original atom and if so assigns that matrix element the
    smallest combination of connectors.

    Parameter
    ---------
    in_mat : 2-d array-like
        Connectivity matrix to be expanded

    Returns
    -------
    out_mat : 2-d array-like
        Once-expanded connectivity matrix

    """
    out_mat = np.copy(in_mat)
    for i, row in enumerate(in_mat):
        # indices of unconnected atoms
        dangles = []
        # indices of connected atoms
        connectors = []
        for j, element in enumerate(row):
            if element == 0 and j > i:
                dangles.append(j)
            elif element != 0:
                connectors.append(j)
        for dangle in dangles:
            orders = []
            for k, dangle_element in enumerate(in_mat[dangle]):
                if dangle_element != 0 and k in connectors:
                    orders.append(dangle_element + in_mat[i][k])
                # if orders is not empty
                if orders:
                    out_mat[i][dangle] = min(orders)
                    out_mat[dangle][i] = min(orders)
    return out_mat


def complete_expand(in_mat):
    """Expand a matrix until it stops expanding."""
    mat = np.copy(in_mat)
    i = 1
    while True:
        i += 1
        temp_mat = expand_connect(mat)
        if np.array_equal(mat, temp_mat):
            break
        mat = np.copy(temp_mat)

    return mat

def get_connectivity_mat(in_mol):
    """Return the connectivity matrix of the Mol"""
    first_connect = detect_1_connect(in_mol)
    connect_mat = complete_expand(first_connect)

    return connect_mat

def charged_kinds(in_atoms, in_kinds):
    """
    Get charged atom kinds from charged atoms and kinds.

    For each kind of atom to be charged, goes through the list of atoms and
    makes an average of the partial atomic charge of atoms of that type.

    Parameters
    ----------
    in_atoms : Mol object
        The atoms should be charged and some of them at least should be of the
        relevant kind
    in_kinds : list of tuples
        The tuples are of the form (a,b) where a is an element string (like 'C')
        and b is a frozenset of ((element string,order of connection),amount of
        connections). (a,b) is known as an atom kind

    Returns
    -------
    q_kinds : list of tuples
        Each tuple is now (average charge,kind). This tuple is known as a
        charged kind

    """
    q_kinds = []
    for kind in in_kinds:
        charges = []
        for atom in in_atoms:
            if atom.kind == kind:
                charges.append(atom.q)
        if charges:  # if not empty
            avg_charge = sum(charges) / float(len(charges))
        else:
            avg_charge = 0
        q_kinds.append((avg_charge, kind))
    return q_kinds


def assign_charges(char_atoms, unchar_atoms):
    """
    Assign charges from one list of atoms to another list of atoms.

    This is based on the connectivity of the charged atoms, as defined by a
    maximum bond length. The function works for periodic or non periodic systems
    in the input and the output. The uncharged atoms are changed and there is no
    output.

    Parameters
    ----------
    char_atoms : Mol object
        Atoms which already have assigned charge
    char_vectors : 3x3 array-like or None
        Lattice vectors of the input system if it is periodic. If not use None
    unchar_atoms : Mol object
        Atoms which need charges assigned to them
    unchar_vectors : 3x3 array-like or None
        See char_vectors

    """
    # detect the charged atom's connectivity matrix
    char_first = detect_1_connect(char_atoms)
    char_cnct = complete_expand(char_first)

    # get charged atom kinds as a result
    kinds = []
    for i, atom in enumerate(char_atoms):
        atom.set_connectivity(char_atoms, char_cnct[i])
        kinds.append(atom.kind)
    kinds = set(kinds)
    q_kinds = charged_kinds(char_atoms, kinds)

    # detect uncharged atom connectivity
    unchar_first = detect_1_connect(unchar_atoms)
    unchar_cnct = complete_expand(unchar_first)

    # determine kind and cross check with charged kinds
    for i, atom in enumerate(unchar_atoms):
        atom.set_connectivity(unchar_atoms, unchar_cnct[i])
        for q_kind in q_kinds:
            if atom.kind == q_kind[1]:
                atom.q = q_kind[0]
    return


def main(in_xyz, in_log, target, output, bonding, thresh, kind):
    if(in_xyz):
        mol = rf.mol_from_file(in_xyz)
    else:
        mol = rf.mol_from_gauss(in_log)
    charges = rf.read_g_char(in_log, kind)[0]
    cluster = rf.mol_from_file(target)

    mol.set_bonding(bonding=bonding, thresh=thresh)
    cluster.set_bonding(bonding=bonding, thresh=thresh)

    for atom, char in zip(mol, charges):
        atom.q = char

    assign_charges(mol, cluster)

    # warning if some atoms have not been assigned or if some original charges
    # were 0
    bad_atoms = []
    for atom in cluster:
        if abs(atom.q) <= 0.000:
            bad_atoms.append(atom)
    if len(bad_atoms) > 0:
        print("WARNING: " + str(len(bad_atoms)) + " atoms have null charge!")
        print(bad_atoms)

    out_file = open(output, "w")
    out_file.write(str(len(cluster)) + "\n\n")
    for atom in cluster:
        out_file.write(str(atom) + "\n")
    out_file.close()

if __name__ == '__main__':
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("in_log", help="Input .log file with RESP analysis",
                        default="gaussian.log")
    parser.add_argument("target", help="Target .xyz file to assign charges to",
                        default="cluster.xyz")
    parser.add_argument(
        "-i", "--in_xyz", help="Input .xyz file of single molecule if the geometry in the log file is not good")
    parser.add_argument("-o", "--output", help="Name of the output file",
                        default="out_char", type=str)
    parser.add_argument("-b", "--bonding", help="Bonding type to be evaluated. Use 'dis' (default), 'cov' or 'vdw' to start calculating the distance at the centre of the atoms, the surface of the covalent sphere or the surface of the vdw sphere.", default="dis", type=str)
    parser.add_argument("-t", "--threshold", help="Maximum length in Angstrom that qualifies as a bond. Default 1.7",
                        default=1.7, type=float)
    parser.add_argument("-k", "--kind", help="Kind of population, mulliken or esp",
                        default="esp", type=str)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    main(args.in_xyz, args.in_log, args.target,
         args.output, args.bonding, args.threshold, args.kind)
