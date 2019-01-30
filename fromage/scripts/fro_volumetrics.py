#!/usr/bin/env python3
import argparse
import fromage.utils.volume as vo
import fromage.io.read_file as rf
import numpy as np
import time
import sys

def main():
    start = time.time()
    in_atoms = sys.argv[1]


    out_file = open("volumes","w")
    prox_grid = vo.CubeGrid()
    vdw_grid = vo.CubeGrid()

    atoms = rf.mol_from_file(in_atoms)
    mol = atoms.select(int(sys.argv[2])-1)
    c_x, c_y, c_z = mol.centroid()
    prox_grid.grid_from_point(c_x, c_y, c_z, res=100,box=np.array([[25.0, 0.0, 0.0], [0.0,25.0,0.0], [0.0, 0.0, 25.0]]))
    vdw_grid.grid_from_point(c_x, c_y, c_z, res=100,box=np.array([[25.0, 0.0, 0.0],[0.0, 25.0,0.0], [0.0, 0.0, 25.0]]))

    rest = []
    for atom in atoms:
        if atom not in mol:
            rest.append(atom)

    prox_grid.set_grid_coord()
    vdw_grid.set_grid_coord()

    prox_grid.proximity(mol,rest)
    vdw_grid.vdw_vol(mol)

    prox_grid.out_cube("prox.cube", atoms)
    out_file.write("Proximity volume: "+str(prox_grid.volume())+"\n")

    vdw_grid.out_cube("vdw.cube", atoms)
    out_file.write("VDW volume: "+str(vdw_grid.volume())+"\n")

    prox_grid.add_grid(vdw_grid.grid)
    out_file.write("Union volume: "+str(prox_grid.volume())+"\n")
    prox_grid.out_cube("add.cube", atoms)

    end = time.time()
    out_file.write("Time Elapsed: "+ str(end - start)+"s")

    out_file.close()



if __name__ == '__main__':
    start = time.time()
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "in_xyz", help="Input .xyz file with the unit cell geometry")
    parser.add_argument(
        "atom_label", help="Label of atom in the molecule to analyse", default=1)
    parser.add_argument(
        "-o", "--output", help="Name of the output file", default="out_cell.xyz", type=str)
    parser.add_argument(
        "-b", "--bonding", help="Type of bonding used for detecting full molecules. The options are dis, cov and vdw", default="dis", type=str)
    parser.add_argument("-T", "--thresh", help="Threshold distance to select a bond. The default pairs are dis:1.8, cov:0.2, vdw:-0.3. To select default, just use a threshold of 999", default=999, type=float)
    parser.add_argument(
        "-c", "--complete", help="Complete the molecules in the cell. Incompatible with -C or -f", action="store_true")
    parser.add_argument(
        "-C", "--confine", help="Confine the molecule to the primitive cell. Incompatible with -c or -f", action="store_true")
    parser.add_argument(
        "-f", "--fractional", help="Print the xyz in fractional coordinates. Incompatible with -c or -C", action="store_true")
    parser.add_argument(
        "-m", "--mono", help="Boolean to print all unique monomers, requires -c. NB: This breaks severely if the wrong -l is chosen", action="store_true")
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
         args.remove_duplicate_atoms, args.output, args.bonding, args.thresh, args.mono, args.translations, args.radius, args.inclusivity, args.center)
    end = time.time()
    print("\nTotal time: {}s".format(round((end - start), 1)))
