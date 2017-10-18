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

in_f = sys.argv[1]

atoms = rf.read_xyz(in_f)[-1]
natoms= len(atoms)
selected=[]

for i in range(natoms):
    mol=ha.select(1.7,atoms,i)
    if mol[0] in [val for sublist in selected for val in sublist]:
        next
        #raise ValueError("Atom " + str(i) + " was already selected!")

    else:
        selected.append(mol)

for i in selected:
    print "{}\n".format(len(i))
    for j in i:
        print j
