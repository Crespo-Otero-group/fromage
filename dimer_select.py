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
    difference=np.array(A) - np.array(B)
    SAD=np.sum(np.abs(np.array(A) - np.array(B)))
    SSD=np.sum(np.square(np.array(A) - np.array(B)))
    correlation=np.corrcoef(np.array((A, B)))[0, 1]
    """print "difference:", difference
    print "SAD:", SAD
    print "SSD:", SSD
    print "correlation:",correlation"""
    return float(SSD)

def bond_distances(dimers):
    connect_mats=[]
    for dim_no,dim_atom in enumerate(dimers):
        dim_cons=[]
        for atom_no_A,atom_A in enumerate(dim_atom):
            for atom_no_B,atom_B in enumerate(dim_atom[atom_no_A:]):
                if atom_A!=atom_B:
                    dim_cons.append(round(vector_distance((atom_A.x,atom_A.y,atom_A.z,atom_B.x,atom_B.y,atom_B.z)),0))
        connect_mats.append(sorted(dim_cons))
    return connect_mats


in_f = sys.argv[1]

atoms = rf.read_xyz(in_f)[-1]
natoms= len(atoms)
print "{} atoms".format(natoms)

##### SELECT MOLECULE
print "Generating molecules..."
selected=make_molecules(atoms,1.7)
print "{} molecules generated".format(len(selected))

###### SELECT DIMERS
print "Generating dimers..."
dimers=make_dimers(selected)
print "{} dimers generated".format(len(dimers))

####### SELECT UNIQUE DIMERS

print "Finding unique dimers..."
distances=bond_distances(dimers)
if len(distances)==1: #dangerous!! Check it works
    ef.write_xyz(sys.argv[1][:-4]+"_unique.xyz",dimers[0])
    exit("One unique dimer found, writing to xyz")
else:
    unique_dims=[]
    unique_distances=[]
    for i,j in enumerate(distances):
        for k,l in enumerate(distances):
            if differences(j,l)>1 and (j not in unique_distances and l not in unique_distances):
                unique_distances.append(j)
                unique_distances.append(l)
                unique_dims.append(dimers[i])
                unique_dims.append(dimers[k])
                print "{} {} : {}".format(i,k,differences(j,l))

print "Number of unique dimers: {}\n".format(len(unique_dims))

for dim_no,dim in enumerate(unique_dims):
        ef.write_xyz(str(sys.argv[1][:-4]+"_unique_"+str(dim_no)+".xyz"),dim)
alldims = [item for sublist in unique_dims for item in sublist]
ef.write_xyz("alldims.xyz",alldims)
end = time.time()
print "Total time: {}s".format(round((end - start),3))
