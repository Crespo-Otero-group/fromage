#!/usr/bin/env python
"""
Utility for selecting unique dimers from a .xyz file

The unique dimers are written to separate output files, *_dimer_*.xyz

The routines herein are built upon implementations of the dimer_select.py module
by Michael Dommett (see earlier commits) and collected personal scripts by
Amir Sidat.

"""

from __future__ import division
import sys
import argparse
import time
import fromage as fro
import numpy as np
from fromage.io import read_file as rf

# safe printing
def prints(arg_in):
    if args.verbose:
        print(arg_in)
    return

def all_dimers(mono_list):
    """
    Generate all possible dimers from a list of Mols
    """
    out_dims = []
    for i, mono_i in enumerate(mono_list[:-1]):
        for mono_j in mono_list[i+1:]:
            dim = fro.make_dimer(mono_i,mono_j)
            out_dims.append(dim)
    return  out_dims

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input .xyz file", type=str)
    parser.add_argument(
        "-vec","--vectors", help="Input lattice vectors (optional)", default="")
    parser.add_argument(
        "-b", "--bonding", help="Type of bonding used for detecting full molecules. The options are dis, cov and vdw", default="cov", type=str)
    parser.add_argument(
        "-t", "--thresh", help="Threshold distance to select a bond. The default pairs are dis:1.8, cov:0.2, vdw:-0.3", default=None, type=float)
    parser.add_argument(
        "-bs", "--bonding_string", help="Alternate specification of bonding. Here the threshold and bonding are lumped up in one string like 'cov-0.1' or 12dis'", default="", type=str)
    parser.add_argument("-T", "--dimtype", help="Use centroid distance 'centroid' or shortest atomic distances 'dis' or covalent radius 'cov' or van der waals radii 'vdw' to define a dimer",
                        default=str("centroid"), type=str)
    parser.add_argument("-d", "--dist", help="Distance criterion (in units of Angstrom) to define a dimer",
                        default=7, type=float)
    parser.add_argument(
        "-lin", "--linear", help="The molecule is linear-like. This means that the principal axis will simply be the longest interatomic axis.", action="store_true")

    parser.add_argument(
        "-nh", "--no_hydrogens", help="Ignore hydrogens in the selection and analysis of dimers", action="store_true")
    parser.add_argument(
        "-na", "--no_atom_label", help="Ignore atoms of the same kind as the selected ones in the selection and analysis of dimers", default=[], type=int, nargs='*')
    parser.add_argument(
        "-l", "--tol_duplicate", help="RMSD tolerance for two dimers to be considered the same. Default 0.0001", default=10e-4, type=float)
    parser.add_argument(
        "-p", "--print_dimers", help="Write out each dimer in .xyz format", action="store_true")
    parser.add_argument(
        "-o", "--output_geometry_data", help="Output file for geometry data analysis", default='dimers.dat', type=str)
    parser.add_argument(
        "-v", "--verbose", help="Print full output", action="store_true")

    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    return args

def test():
    from fromage.utils.atom import Atom
    from fromage.utils.mol import Mol
    mol_a = Mol([Atom("C",1)])
    mol_b = Mol([Atom("C",2)])
    mol_c = Mol([Atom("C",3)])
    mol_d = Mol([Atom("C",4)])

    lis = [mol_a,mol_b,mol_c,mol_d]
    print(all_dimers(lis))

def main(args):
    # read input file
    all_atoms = fro.mol_from_file(args.input)
    # specify bonding
    all_atoms.set_bonding(bonding=args.bonding, thresh=args.thresh)
    # if the bonding is specified all in one string
    if args.bonding_string:
        all_atoms.set_bonding_str(args.bonding_string)

    # if periodic
    if args.vectors:
        prints("Vectors detected")
        vectors = rf.read_vectors(args.vectors)
        all_atoms.vectors = vectors
        # if we are forbidding some atoms
        if args.no_atom_label:
            forbidden_kinds = []
            # for each atom to forbid
            for label in args.no_atom_label:
                mol_of_interest = all_atoms.per_select(label-1)
                mol_of_interest.set_connectivity()
                for atom in mol_of_interest:
                    if atom == all_atoms[label-1]:
                        forbidden_kinds.append(atom.kind)
        all_atoms, molecules = all_atoms.complete_cell()
        prints("{} molecules detected after reconstitution".format(len(molecules)))
    else:
        # if we are forbidding some atoms
        if args.no_atom_label:
            forbidden_kinds = []
            # for each atom to forbid
            for label in args.no_atom_label:
                mol_of_interest = all_atoms.select(label-1)
                mol_of_interest.set_connectivity()
                for atom in mol_of_interest:
                    if atom == all_atoms[label-1]:
                        forbidden_kinds.append(atom.kind)
        # separate into molecules
        molecules = all_atoms.segregate(diff_mols=False)
        prints("{} molecules detected".format(len(molecules)))

    # ignore hydrogens
    if args.no_hydrogens:
        for mol in molecules:
            mol.geom.ignore_hydrogens = True

    # ignore certain atoms
    if args.no_atom_label:
        zero_vectors = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
        for mol in molecules:
            mol.vectors = zero_vectors
            mol.set_connectivity()
            mol.geom.ignore_kinds = forbidden_kinds

    # get dimers
    dimers_raw = all_dimers(molecules)
    prints("{} dimers in the input geometry".format(len(dimers_raw)))

    # if periodic
    if args.vectors:
        # get all the dimers generated by extra periodic images
        dimers_per = []
        for dimer in dimers_raw:
            images = dimer.images(vectors)
            dimers_per.extend(images)
        dimers_raw = dimers_per
        prints("{} dimers considering periodicity".format(len(dimers_raw)))

    if args.linear:
        for dimer in dimers_raw:
            dimer.mols_are_linear()
    # filter out dimers
    if args.dimtype == 'centroid':
        method = 'centroid'
        mode = 'dis'
    else:
        method = 'atomic'
    if args.dimtype == 'dis':
        mode = 'dis'
    elif args.dimtype == 'cov':
        mode = 'cov'
    elif args.dimtype == 'vdw':
        mode = 'vdw'

    selected_dimers = []
    for dimer in dimers_raw:
        if dimer.inter_distance(method=method,mode=mode) <= args.dist:
            selected_dimers.append(dimer)
    prints("{} dimers selected from distance criterion".format(len(selected_dimers)))

    # remove repeated dimers
    keep_these = []
    for i,dimer_i in enumerate(selected_dimers[:-1]):
        keep = True
        for dimer_j in selected_dimers[i+1:]:
            if dimer_i.same_geom(dimer_j,tol=args.tol_duplicate):
                keep = False
                break
        if keep:
            keep_these.append(dimer_i)
    # don't forget that the last one was not included in the parent loop
    keep_these.append(selected_dimers[-1])
    selected_dimers = keep_these
    prints("{} dimers remaining after removing duplicates".format(len(selected_dimers)))

    if args.output_geometry_data:
        data_file = open(args.output_geometry_data, "w")
        header_str = "{:>13}{:>10}{:>10}{:>10}{:>17}{:>13}{:>11}\n".format("Dimer number","Alpha","Beta","Gamma","Centroid dist","Slip angle","Category")
        data_file.write(header_str)
    for i,dimer in enumerate(selected_dimers):

        if args.output_geometry_data:
            cen_dist = dimer.inter_distance(method='centroid')
            angles = dimer.angles()
            slip_angle = dimer.slip_angle()
            # classify
            if 20 < angles[2] < 160:
                category = "E-F"
            elif slip_angle > 60:
                category = "S-S"
            else:
                category = "F-F"
            data_file.write("{:>7}{:17.3f}{:10.3f}{:10.3f}{:12.3f}{:15.3f}{:>10}\n".format(i+1,angles[0],angles[1],angles[2],cen_dist,slip_angle,category))
        if args.print_dimers:
            out_name = str(args.input[:-4]) + "_dimer_" + str(i+1) + ".xyz"
            dimer.write_xyz(out_name)
    if args.output_geometry_data:
        data_file.close()

if __name__ == '__main__':
    start = time.time()
    args = parse_args()
    main(args)
    end = time.time()
    prints("\nTotal time: {}s".format(round((end - start), 1)))
