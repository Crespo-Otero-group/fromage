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
    """
    Calculate the distances between two cartesian coordinates

    Parameters
    ----------
    (x1,y1,z1,x2,y2,z2): hextuple of floats
            atomic coordinates of atoms 1 and 2
    Returns
    -------
    dist: Float
            distance in units of coordinates
    """
    #calculate distance
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2-z1)**2)
    return dist

def make_molecules(atoms,bl):
    """
    Generate list of molecules based on bond length bl

    Parameters
    ----------
    atoms: list of Atom objects
    bl: float
        Bond length in unit of input file
    Returns
    -------
    selected: list of lists
        List L of length M molecules, where each member of L is a list of atom objects
    """
    selected=[]
    max_length=0
    for i,atom in enumerate(atoms):
        if atom not in [val for sublist in selected for val in sublist]:
            molecule=ha.select(bl,atoms,i) #creates a molecule
            if len(molecule)>max_length:
                max_length=len(molecule)
                selected=[]
                selected.append(molecule)
            elif len(molecule)==max_length:
                selected.append(molecule)

    return selected

def make_dimers(selected,cd):
    """
    Generate a list of dimers based on centroid distances cd

    Parameters
    ----------
    selected: list of lists
        M molecules containing N atom objects
    cd: float
        Length between centroids of molecules
    Returns
    -------
    dimers: list of lists
        List L of length D dimers, where each member of L is a list of 2N atom objects
    """
    dimers=[]
    for mol_1_no,mol1 in enumerate(selected):
        for mol_2_no,other_mol in enumerate(selected[mol_1_no:]):
            if mol1!=other_mol:
                cent_1=ha.find_centroid(mol1)
                cent_2=ha.find_centroid(other_mol)
                if vector_distance(cent_1+cent_2) <=cd:
                    vector_distance(cent_1+cent_2)
                    new_mol=mol1+other_mol
                    dimers.append(new_mol)
    return dimers

def differences(A,B):
    """
    Calulate the sum of squares difference between two lists, nominally of atomic distances

    Parameters
    ----------
    A,B: lists of floats
        Bond distances in dimers
    Returns
    -------
    SSD: float
        Sum of squares differences
    """
    SSD=np.sum(np.square(np.array(A) - np.array(B)))
    return float(SSD)

def interatomic_distances(dimers):
    """
    Generate a list of the interatomic distances for each dimer in list of dimers

    Parameters
    ----------
    dimers: list of lists
     List of L of length D dimers, where each member of L is a list of 2N atom objects
    Returns
    -------
    connections: list
        Connections per dimer
    """
    connections=[]
    for dim_no,dim_atom in enumerate(dimers):
        dim_cons=[]
        for atom_no_A,atom_A in enumerate(dim_atom):
            for atom_no_B,atom_B in enumerate(dim_atom[atom_no_A:]):
                if atom_A!=atom_B:
                    dim_cons.append(round(vector_distance((atom_A.x,atom_A.y,atom_A.z,atom_B.x,atom_B.y,atom_B.z)),0))
        connections.append(sorted(dim_cons))
    return connections

if __name__ == "__main__":
    # parse the input
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", help="Input .xyz file")
    parser.add_argument("-b", "--bond", help="Maximum length (in unites of input file) that qualifies as a bond",
                        default=1.6, type=float)
    parser.add_argument("-c","--centdist", help="Distance criterion (in unites of input file) to define a dimer, between the cetroids of two monomers",
                        default=7.0, type=float)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    atoms = rf.read_xyz(args.input)[-1]
    natoms= len(atoms)
    print "{} atoms".format(natoms)

    ##### SELECT MOLECULE
    print "Generating molecules..."
    selected=make_molecules(atoms,args.bond)
    print "{} molecules generated".format(len(selected))

    ###### SELECT DIMERS
    print "Generating dimers..."
    dimers=make_dimers(selected,args.centdist)
    print "{} dimers generated".format(len(dimers))

    ####### SELECT UNIQUE DIMERS

    print "Finding unique dimers..."
    distances=interatomic_distances(dimers)
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


    print "Number of unique dimers: {}".format(len(unique_dims))

    # write the files
    for dim_no,dim in enumerate(unique_dims):
            outfile=str(args.input[:-4])+"_dimer_"+str(dim_no)+".xyz"
            print "Writing {}".format(outfile)
            ef.write_xyz(outfile,dim)

    end = time.time()
    print "Total time: {}s".format(round((end - start),3))
