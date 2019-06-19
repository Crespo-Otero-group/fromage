#!/usr/bin/env python3
import argparse
import fromage.utils.volume as vo
import fromage.io.read_file as rf
import numpy as np
import time
import sys

def parse_args():
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_xyz", help="Input .xyz file with the unit cell geometry")
    parser.add_argument(
        "-l","--atom_label", help="Label of atom in the molecule to analyse", default=1, type=int)
    parser.add_argument(
        "-b", "--bonding", help="Type of bonding used for detecting full molecules. The options are dis, cov and vdw", default="dis", type=str)
    parser.add_argument(
        "-T", "--thresh", help="Threshold distance to select a bond. The default pairs are dis:1.8, cov:0.2, vdw:-0.3. To select default, just use a threshold of 999", default=999, type=float)
    parser.add_argument(
        "-bs", "--bonding_string", help="Alternate specification of bonding. Here the threshold and bonding are lumped up in one string like 'cov-0.1' or 12dis'", default="", type=str)
    parser.add_argument("-res", "--resolution",
                        help="The number of voxels per side of the box", default=100, type=int)
    parser.add_argument("-dim", "--dimensions",help="Dimensions of the box in Angstrom. Give x y and z like '30 30 30'",default=[25, 25, 25], type=float, nargs='*')

    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)

    return args
def main(args):

    in_atoms = args.in_xyz
    atoms = rf.mol_from_file(in_atoms)

    out_file = open("volumes", "w")
    prox_grid = vo.CubeGrid()
    vdw_grid = vo.CubeGrid()

    mol = atoms.select(args.atom_label - 1)
    c_x, c_y, c_z = mol.centroid()

    x_dim, y_dim, z_dim = args.dimensions[0], args.dimensions[1], args.dimensions[2]

    prox_grid.grid_from_point(c_x, c_y, c_z, res=args.resolution, box=np.array(
        [[x_dim, 0.0, 0.0], [0.0, y_dim, 0.0], [0.0, 0.0, z_dim]]))
    vdw_grid.grid_from_point(c_x, c_y, c_z, res=100, box=np.array(
        [[x_dim, 0.0, 0.0], [0.0, y_dim, 0.0], [0.0, 0.0, z_dim]]))

    rest = []
    for atom in atoms:
        if atom not in mol:
            rest.append(atom)

    prox_grid.set_grid_coord()
    vdw_grid.set_grid_coord()

    prox_grid.proximity(mol, rest)
    vdw_grid.vdw_vol(mol)

    prox_grid.out_cube("voro.cube", atoms)
    out_file.write("Voronoi volume: " + str(prox_grid.volume()) + "\n")

    vdw_grid.out_cube("vdw.cube", atoms)
    out_file.write("VDW volume: " + str(vdw_grid.volume()) + "\n")

    prox_grid.add_grid(vdw_grid.grid)
    out_file.write("Union volume: " + str(prox_grid.volume()) + "\n")
    prox_grid.out_cube("add.cube", atoms)

    out_file.close()


if __name__ == '__main__':
    start = time.time()
    args = parse_args()
    main(args)
    end = time.time()
    print("\nTotal time: {}s".format(round((end - start), 1)))
