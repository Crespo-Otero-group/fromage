#!/usr/bin/env python
"""
Utility for selecting unique dimers from a .xyz file

The default output is called out.xyz but can be specified with -o

"""
import time

start = time.time()
import sys
import argparse
import assign_charges as ac
import read_file as rf
import edit_file as ef
import handle_atoms as ha
from math import sqrt
import numpy as np

def vector_distance((x1,y1,z1,x2,y2,z2)):
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2-z1)**2)
    return dist

def make_molecules(atoms,bl):
    selected=[]
    max_length=0
    for i,atom in enumerate(atoms):
        if atom not in [val for sublist in selected for val in sublist]:
            molecule=ha.select(1.7,atoms,i) #creates a molecule
            if len(molecule)>max_length:
                max_length=len(molecule)
                selected=[]
                selected.append(molecule)
            elif len(molecule)==max_length:
                selected.append(molecule)

    return selected

def make_dimers(selected):
    dimers=[]
    for mol_1_no,mol1 in enumerate(selected):
        for mol_2_no,other_mol in enumerate(selected[mol_1_no:]):
            if mol1!=other_mol:
                cent_1=ha.find_centroid(mol1)
                cent_2=ha.find_centroid(other_mol)
                if vector_distance(cent_1+cent_2) <7:
                    new_mol=mol1+other_mol
                    dimers.append(new_mol)
    return dimers

def differences(A,B):
    print "difference:", A - B
    print "SAD:", np.sum(np.abs(A - B))
    print "SSD:", np.sum(np.square(A - B))
    print "correlation:", np.corrcoef(np.array((A, B)))[0, 1]

in_f = sys.argv[1]

atoms = rf.read_xyz(in_f)[-1]
natoms= len(atoms)

##### SELECT MOLECULE
print "{} atoms".format(natoms)
print "Making molecules..."
selected=make_molecules(atoms,1.7)
print "{} molecules generated".format(len(selected))

###### SELECT DIMERS

print "Making dimers..."
dimers=make_dimers(selected)
print "{} dimers generated".format(len(dimers))

####### SELECT UNIQUE DIMERS

for i in range(8):
    ef.write_xyz(sys.argv[1][:-4]+"dim-"+str(i)+".xyz",dimers[i])
print "Finding unique dimers..."
connect_mats=[]
for dim_no,dim_atom in enumerate(dimers):
    dim_cons=[]
    for atom_no_A,atom_A in enumerate(dim_atom):
        for atom_no_B,atom_B in enumerate(dim_atom[atom_no_A:]):
            if atom_A!=atom_B:
                dim_cons.append(round(vector_distance((atom_A.x,atom_A.y,atom_A.z,atom_B.x,atom_B.y,atom_B.z)),0))
    connect_mats.append(sorted(dim_cons))
differences(connect_mats[0],connect_mats[1])
different=[]
indexes=[]
if len(connect_mats)==1: #dangerous!! Check it works
    #ef.write_xyz("unique_dimer.xyz",dimers[0])
    exit("One unique dimer found, writing to xyz")
else:
    for i,j in enumerate(connect_mats):
        for k,l in enumerate(connect_mats[i:]):
            if j != l and j not in different:
                different.append(j)
                indexes.append(i)
print indexes

unique_dims=[dimers[i] for i in indexes]
print "Number of unique dimers: {}\n".format(len(unique_dims))

for dim_no,dim in enumerate(unique_dims):
        ef.write_xyz(str(sys.argv[1][:-4]+"_unique_"+str(dim_no)+".xyz"),dim)

alldims = [item for sublist in unique_dims for item in sublist]
ef.write_xyz("alldims.xyz",alldims)
end = time.time()
print "Total time: {}s".format(round((end - start),0))
