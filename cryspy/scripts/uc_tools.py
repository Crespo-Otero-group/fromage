#!/usr/bin/env python
"""
Utility to manipulate unit cell geometries in xzy format to get information like
unique monomer identities and modified 'uncropped' unit cells.

"""
from cryspy.io import read_file as rf
from cryspy.io import edit_file as ef
from cryspy.io import parse_config_file as pcf
from cryspy.utils import handle_atoms as ha

import numpy as np

def main(in_xyz, vectors_file, output, max_r):
    vectors = rf.read_vectors(vectors_file)
    atoms = read_pos(in_xyz)
    mod_cell = ha.complete_cell(atoms,vectors,max_bl=max_r)


print(vectors)

if __name__ == '__main__':
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("in_xyz", help="Input .xyz file with the unit cell geometry")
    parser.add_argument("vectors", help="Input lattice vectors", default="vectors")
    parser.add_argument("-o", "--output", help="Name of the output file",default="out_cell.xyz", type=str)
    parser.add_argument("-b", "--bond", help="Maximum length in Angstrom that qualifies as a bond. Default 1.7",
                        default=1.7, type=float)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    main(args.in_xyz, args.vectors, args.output, args.bond)
