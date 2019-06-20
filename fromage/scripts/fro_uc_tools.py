#!/usr/bin/env python
"""
Utility to manipulate unit cell geometries in xzy format to get information like
unique monomer identities and modified 'uncropped' unit cells.

"""
import argparse
import numpy as np
import sys
from copy import deepcopy
import time

from fromage.io import read_file as rf
from fromage.io import edit_file as ef


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


def main(in_xyz, vectors_file, complete, confine, frac, dupli, output, bonding, thresh, bonding_str, print_mono, trans, clust_rad, inclusivity, center_label):
    vectors = rf.read_vectors(vectors_file)
    atoms = rf.mol_from_file(in_xyz,vectors=vectors)

    atoms.set_bonding(bonding=bonding, thresh=thresh)
    # if the bonding is specified all in one string
    if bonding_str:
        atoms.set_bonding_str(bonding_str)

    if print_mono and not complete:
        print("-c is required for -m")
        return

    if sum([complete, confine, frac]) > 1:
        print("-c, -C and -f are mutually exclusive")
        return

    print_now = True
    centered = False
    if center_label:
        central_mol, atoms = atoms.centered_mols(center_label-1) # -1 for Python index
        centered = True
    if complete:
        # get the modified (uncropped) unit cell
        mod_cell, mols = atoms.complete_cell()
    elif confine:
        # get the cell confined to primitive
        mod_cell = atoms.confined()
    elif frac:
        # get the cell in fractional coordinates
        mod_cell = atoms.dir_to_frac_pos()
    elif centered:
        mod_cell = atoms.copy()
    else:
        mod_cell = deepcopy(atoms)
        print_now = False

    if dupli:
        # purge duplicate atoms
        mod_cell.remove_duplicates()
        print_now = True

    if print_now:
        mod_cell.write_xyz(output)

    if trans:
        translation = np.array(trans)
        new_atoms = mod_cell.supercell(translation)
        new_vec = (vectors.T * translation.transpose()).T
        new_atoms.write_xyz("supercell_out.xyz")
        ef.write_lat_vec("supercell_vectors", new_vec)

    if clust_rad > 0.0:
        clust = atoms.make_cluster(clust_rad, mode = inclusivity)
        clust.write_xyz("cluster_out.xyz")

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
                for at_i in mol_i:
                    for at_j in mol_j:
                        dimer_distances.append(at_i.per_dist(at_j,vectors))
                        dimer_distances.sort()
                identity.append(dimer_distances)
                identity.sort()
            identities.append(identity)
        identities = np.array(identities)

        # list of each molecules paired with its identity
        # in Python 2 we could just use a zip but for compatibility with Python3
        # we need a list of tuples
        mols_ids = [(mol_i,iden_i) for mol_i, iden_i in zip(mols, identities)]

        # list of each molecule paired with its identity separated by equivalent
        # molecule
        sep_mols_ids = []
        while len(mols_ids) > 0:
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
            ef.write_xyz("mono_" + str(counter + 1) + ".xyz", fequimols)
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
    parser.add_argument(
        "-b", "--bonding", help="Type of bonding used for detecting full molecules. The options are dis, cov and vdw", default="dis", type=str)
    parser.add_argument(
        "-T", "--thresh", help="Threshold distance to select a bond. The default pairs are dis:1.8, cov:0.2, vdw:-0.3", default=None, type=float)
    parser.add_argument(
        "-bs", "--bonding_string", help="Alternate specification of bonding. Here the threshold and bonding are lumped up in one string like 'cov-0.1' or 12dis'", default="", type=str)
    parser.add_argument(
        "-c", "--complete", help="Complete the molecules in the cell. Incompatible with -C or -f", action="store_true")
    parser.add_argument(
        "-C", "--confine", help="Confine the molecule to the primitive cell. Incompatible with -c or -f", action="store_true")
    parser.add_argument(
        "-f", "--fractional", help="Print the xyz in fractional coordinates. Incompatible with -c or -C", action="store_true")
    parser.add_argument(
        "-m", "--mono", help="Boolean to print all unique monomers, requires -c. NB: This breaks severely if the wrong -bs is chosen", action="store_true")
    parser.add_argument("-t", "--translations",
                        help="Create a supercell via lattice translations", default=None, type=int, nargs='*')
    parser.add_argument("-d", "--remove_duplicate_atoms",
                        help="Purge duplicate atoms", action="store_true")
    parser.add_argument("-r", "--radius", help="Generate a cluster of molecules of the given radius. Radius 0.0 turns this off.",
                        default=0.0, type=float)
    parser.add_argument("-i", "--inclusivity", help="Choose between inclusive (inc) or exclusive (exc) radius selecting.", default='exc', type=str),
    parser.add_argument("-e", "--center", help="Move the atoms of the cell such that the molecule containing the specified label is at the origin. To turn off: 0",default=0,type=int)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    main(args.in_xyz, args.vectors, args.complete, args.confine, args.fractional,
         args.remove_duplicate_atoms, args.output, args.bonding, args.thresh, args.bonding_string, args.mono, args.translations, args.radius, args.inclusivity, args.center)
    end = time.time()
    print("\nTotal time: {}s".format(round((end - start), 1)))
