#!/usr/bin/env python
######################################################
#
# fromage module for scaling wigner sampled structure
#
# Reproduced with the permision of PYRAI2MD
#
######################################################
import sys
import numpy as np

def read_xyz(xyz):
    natom = int(xyz[0])
    atoms = [x.split()[0] for x in xyz[2: 2 + natom]]
    coord = np.array([x.split()[1: 4] for x in xyz[2: 2 + natom]]).astype(float)
    return atoms, coord

def write_xyz(atoms, coord, idx):
    natom = len(atoms)
    output = '%s\nGeom %s\n' % (natom, idx)
    for n in range(natom):
        a = atoms[n]
        x, y, z = coord[n]
        output += '%-5s%24.16f%24.16f%24.16f\n' % (a, x, y, z)

    return output

def read_init_xyz(initxyz):
    coord = []
    for n, i in enumerate(initxyz):
        if 'Initial Conditions' in i:
            continue
        if 'NATOM' in i:
            natom = int(i.split()[1])
        if 'Condition' in i:
            initcond = np.array([i.split()[1:] for i in initxyz[n + 1: n + 1 + natom]]).astype(float)
            coord.append(initcond)

    return coord

def write_init_xyz(atoms, coord, idx):
    natom = len(atoms)
    output = 'Condition %5d \n' % (idx)
    for n in range(natom):
        a = atoms[n]
        x, y, z, m, vx, vy, vz = coord[n]
        output += '%-5s%30s%30s%30s%16s%30s%30s%30s\n' % (a, x, y, z, m, vx, vy, vz)

    return output

def main():
    usage = """
        PyRAI2MD initial condition scaling tool

        Usage:
            python3 scaled_wigner.py $nflex_atoms $fac $file.xyz $file.init.xyz

        """

    if len(sys.argv) <= 3:
        exit(usage)

    nflex_atoms = int(sys.argv[1])
    title = sys.argv[2].split('.')[0]
    fac = float(sys.argv[2])
    with open(sys.argv[3], 'r') as infile:
        ref = infile.read().splitlines()

    with open(sys.argv[4], 'r') as infile:
        initxyz = infile.read().splitlines()

    atoms, ref = read_xyz(ref)
    initxyz = read_init_xyz(initxyz)

    out_xyz = ''
    out_initxyz = ''
    for n, xyz in enumerate(initxyz):
        xyz[:nflex_atoms, 0: 3] = (xyz[:nflex_atoms, 0: 3] - ref[:nflex_atoms]) * fac + ref[:nflex_atoms]
        xyz[:, 4: 7] = xyz[:, 4: 7] * fac
        out_xyz += write_xyz(atoms, xyz[:, 0: 3], n + 1)
        out_initxyz += write_init_xyz(atoms, xyz, n + 1)

    with open('%s-scaled-%s.xyz' % (title, fac), 'w') as out:
        out.write(out_xyz)

    with open('%s-scaled-%s.init.xyz' % (title, fac), 'w') as out:
        out.write(out_initxyz)


if __name__ == '__main__':
    main()
