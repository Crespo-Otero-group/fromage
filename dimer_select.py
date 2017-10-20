#!/usr/bin/env python
"""
Utility for selecting unique dimers from a .xyz file

The default output is called out.xyz but can be specified with -o

"""
import sys
import argparse
import read_file as rf
import edit_file as ef
import handle_atoms as ha
from math import sqrt

def vector_distance((x1,y1,z1,x2,y2,z2)):
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2-z1)**2)
    return dist

in_f = sys.argv[1]

atoms = rf.read_xyz(in_f)[-1]
natoms= len(atoms)
selected=[]


for i,atom in enumerate(atoms):
    if atom in [val for sublist in selected for val in sublist]:
        next # skip atom if already assigned
    else:
        molecule=ha.select(1.7,atoms,i) #creates a molecule
        selected.append(molecule)

dimers=[]
for mol1 in selected:
    cent_1=ha.find_centroid(mol1)
    for other_mol in selected:
        if other_mol != mol1:
            cent_2=ha.find_centroid(other_mol)
            if vector_distance(cent_1+cent_2) <6.55:
                if other_mol in dimers:
                    next
                else:
                    print vector_distance(cent_1+cent_2
                    new_mol=mol1+other_mol
                    dimers.append(mol1+other_mol)

        else:
            next
print len(dimers)
