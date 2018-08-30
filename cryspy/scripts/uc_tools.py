#!/usr/bin/env python
"""
Utility to manipulate unit cell geometries in xzy format to get information like
unique monomer identities and modified 'uncropped' unit cells.

"""
import argparse
import numpy as np
import sys
from copy import copy
import itertools as it
import time

from cryspy.io import read_file as rf
from cryspy.io import edit_file as ef
from cryspy.utils import handle_atoms as ha


def compare_id(id_i, id_j):
    """
    Compare two sets of dimeric distance identities

    Parameters
    ----------
    id_i, id_j : list of list of floats
        Each identity contains a list of dimeric distances for each dimer involving
        the relevant molecule in the system. Each dimeric distances is a list of
        intermolecular interatomic distances across the dimer.
    Returns
    -------
    sdiff : float or None
        The comparison index

    """
    # If the sorted lists have different dimensions, the identities are probably
    # different
    if all(len(i) != len(j) for i, j in zip(id_j, id_i)):
        sdiff = None
    else:
        sdiffs = []
        denom = 0
        for i, j in zip(id_i, id_j):
            sdiffs.append(np.sum(np.square(i - j)))
            denom += len(i)
        sdiff = np.sum(sdiffs) / denom
#        sdiff = max(sdiffs)/len(id_i[0])
        return sdiff


def main(in_xyz, vectors_file, complete, output, max_r, print_mono, trans):
    vectors = rf.read_vectors(vectors_file)
    atoms = rf.read_pos(in_xyz)

    if print_mono and not complete:
        print("-c is required for -m")
        return
    if complete:
        # get the modified (uncropped) unit cell
        mod_cell, mols = ha.complete_cell(atoms, vectors, max_bl=max_r)
        ef.write_xyz(output, mod_cell)
    else:
        mod_cell = copy(atoms)

    if trans:
        new_atoms,new_vec = ha.supercell(mod_cell, vectors, np.array(trans))
        ef.write_xyz("supercell_out.xyz",new_atoms)
        ef.write_lat_vec("supercell_vectors",new_vec)

    if print_mono:
        identities = []
        # for each molecule of the modified unit cell
        for mol_i in mols:
            # get the rest of the molecules
            rest = []
            for mol_j in mols:
                if mol_j != mol_i:
                    rest.append(mol_j)
            identity = []
            # for each dimer including the selected molecule
            for mol_j in rest:
                dimer_distances = []
                # get all of the interatomic distances across the dimer using the
                # shortest lattice distance
                for pair in it.product(mol_i, mol_j):
                    dimer_distances.append(pair[0].dist_lat(pair[1].x, pair[1].y, pair[
                                           1].z, vectors[0], vectors[1], vectors[2], order=2)[0])
                    dimer_distances.sort()
                identity.append(dimer_distances)
                identity.sort()
            identities.append(identity)
        identities = np.array(identities)

        # list of each molecules paired with its identity
        mols_ids = zip(mols, identities)

        # list of each molecule paired with its identity separated by equivalent
        # molecule
        sep_mols_ids = []
        while len(mols_ids)>0:
            sep_mol_id = [mols_ids[0]]
            for mol_i, id_i in mols_ids[1:]:
                comp = compare_id(mols_ids[0][1], id_i)
                # if the identities being compared don't have the same dimension,
                # the comparison will return None. This implies that the molecules
                # are not equivalent. This case only happens in crystals with
                # different monomers
                if comp:
                    if comp < 0.001:
                        sep_mol_id.append((mol_i, id_i))
                        mols_ids.remove((mol_i, id_i))
            del mols_ids[0]
            sep_mols_ids.append(sep_mol_id)

        for counter, equimols in enumerate(sep_mols_ids):
            # only the list of lists of atoms i.e. list of molecules which are
            # equivalent
            fequimols = [i[0] for i in equimols]
            # make it a flat list
            fequimols = [inner for outer in fequimols for inner in outer]
            ef.write_xyz("mono_"+str(counter+1)+".xyz",fequimols)
    return


if __name__ == '__main__':
    start = time.time()
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_xyz", help="Input .xyz file with the unit cell geometry")
    parser.add_argument(
        "vectors", help="Input lattice vectors", default="vectors")
    parser.add_argument(
        "-o", "--output", help="Name of the output file", default="out_cell.xyz", type=str)
    parser.add_argument("-b", "--bond", help="Maximum length in Angstrom that qualifies as a bond. Default 1.7",
                        default=1.7, type=float)
    parser.add_argument("-c", "--complete", help="Complete the molecules in the cell", action="store_true")
    parser.add_argument("-m", "--mono", help="Boolean to print all unique monomers, requires -c", action="store_true")
    parser.add_argument("-t", "--translations",help="Create a supercell via lattice translations", default=None, type=int, nargs='*')
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    main(args.in_xyz, args.vectors, args.complete, args.output, args.bond, args.mono, args.translations)
    end = time.time()
    print("\nTotal time: {}s".format(round((end - start), 1)))
