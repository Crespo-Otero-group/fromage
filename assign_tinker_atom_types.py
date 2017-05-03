#!/usr/bin/env python

from sys import argv,exit
from collections import defaultdict
import numpy as np

""" This script assigns forcefield atom types to a testfile,
based on a template file. The template file must have atom
types already assigned, and be in the tinker xyz format.
The tinker format allows for connectivities to be used.
The script reads the template and assigns each atom a vector, 
based on the atom's connectivities. The vector is made up of a 
count for each element in alphabetical order. 
The first vector also includes the base element (the one from which you
are tracking connectivites). 
Each vector is added to by connectivities which are 1,2,3 etc bonds away,until 
each atom type can be uniquely described by a vector.

The testfile is then read and vectors are created for each atom.
When a vector in the testfile matches the template, the testfile 
atom is assigned the corresponding template atom type.

Example : Methane, CH2Cl2

Base elements : C(x1), Cl (x2), H(x2)
Vector = [number of Cs, number of Cls, number of Hs]

Initial vector for C = [122]
Initial vector for each Cl = [110]
Initial vector for each H = [101] 


Remember that for the first vector, the base element is also included.

Often, higher order vectors must be generated, where atoms further away
must be included. The program will loop through connections until a unique 
description for atoms of different type can be identified. 

Important to note: two atoms which are indistinguishable by connectivity but
have different atom types will, of course, not be impossible to assign. 

Enjoy!
Michael Dommett 
m.dommett@qmul.ac.uk
June 2016
""" 




###
# Usage: assign_tinker_atom_types.py template.xyz testfile.xyz outfile.xyz
###
print "\n--- Starting Job ---\n"

# 1) import template, read, save, get number of atoms

template = open(argv[1],"r").read().splitlines()
split_template = (template[0]).split()
natoms_template = int(split_template[0])


# 2) import testfile, read, save, get number of atoms

testfile = open(argv[2],"r").read().splitlines()
split_testfile = (testfile[0]).split()
testfile_natoms = int(split_testfile[0])

# 3) assign outfile
outfile = open(argv[3],"w")

### Functions ###

# count the number of unique atoms (eg C,H,O) in the template 

def count_unique(template,atoms):
    print "Counting unique atoms...\n"
    for line in template[1:]:
        split_line = line.split()
        atom = str(split_line[1])
        atoms.append(atom)
    atoms = np.unique(atoms)
    print "Unique atoms are: \n",
    for i in atoms:
        print "{0} ".format(str(i)),
    print "\n"    
        
    return

# Generate the connectivity of the atom(s) one bond away

def gen_1_connections(file_type,atoms,output_numbs,output_symbs,vector_dict):

    for line in file_type[1:]:

        columns = line.split()
        col_len = len(columns)
        atom_no = int(columns[0])
        atom_type1 = str(columns[1])
        output_numbs[atom_no].append(atom_no)
        output_symbs[atom_no].append(atom_type1)
        for con_a in columns[6:col_len]: 
            con_a = int(con_a)
            output_numbs[atom_no].append(con_a)
            con_line = (int(con_a))
            con_line_split = file_type[con_line].split()
            con_symb = str(con_line_split[1])
            output_symbs[atom_no].append(con_symb)
  
    return
# Create a vector for each atom based on the atoms to which it is connected
def vector_dictionary(output_symbs,vector_dict):
    
    for atom in atoms: 

        atom = str(atom)
        count = 0
        for atom_no, atmsymb in output_symbs.items():
            atom_no = int(atom_no)
            for symb in atmsymb:
                symb = str(symb)
                if atom == symb:
                    count += 1

            vector_dict[atom_no].append(count)
            count = 0

    return

# Generate connectivites for atoms more than one bond away

def gen_n_connections(file_type,atoms,input_numbs,input_symbs,output_numbs,output_symbs,vector_dict):

    for atom_no, con in input_numbs.items():
        for con_a in con:
            con_line = int(con_a)
            con_line_split = file_type[con_line].split()
            con_len = len(con_line_split)
            for con_aa in con_line_split[6:con_len]:
                con_aa = int(con_aa)
                output_numbs[atom_no].append(con_aa)
                con_line = (int(con_aa))
                con_line_split = file_type[con_line].split()
                con_symb = str(con_line_split[1])
                output_symbs[atom_no].append(con_symb)

    return

# Check if the vector created by the vector function can uniquely describe different atom types
def check_vector(template,vector_dict):
    a = []
    equiv_atoms = set()
    for atom_no,vector in vector_dict.items():
        for atmno,vctr in vector_dict.items():
            if vector == vctr and atom_no != atmno:
                atom1_line = template[int(atom_no)].split()
                atom1_type = atom1_line[5]
                atom2_line = template[int(atmno)].split()
                atom2_type = atom2_line[5]
                if atom1_type != atom2_type:
                    a.append(atom_no)
                    equiv_atoms.update(a)

    if equiv_atoms:
        return True
    else:
        return False

# Create a dictionary for the template, which maps vectors to atom types,
# and for the testfile, which maps atom numbers to vectors 
    
def generate_dictionaries(file_type,vector_dict,atom_type_dictionary,atom_no_dictionary):

    vect = {}
    for atom_no, vector in vector_dict.items():
        vect[atom_no]= ''.join(map(str, vector))
    for atom_no,vector in vect.items():
        atom_type_line = file_type[atom_no].split()
        atom_type = atom_type_line[5]
        if file_type == template:
            atom_type_dictionary[vector] = atom_type
        if file_type == testfile:
            atom_no_dictionary[atom_no] = vector

    return

# based on the dictionaries, a mapping from vector to atom number can be made

def print_output(atom_no_dictionary,atom_type_dictionary,testfile):

    outfile.write("{0}  Edited Tinker File from {1} and {2}\n".format(testfile_natoms,argv[1],argv[2]))
    for atom_no,vect1 in atom_no_dictionary.items():
        for vect2,atom_type in atom_type_dictionary.items():
            if vect1 == vect2:
                col = testfile[atom_no].split()
                outfile.write( "{0:>4} {1:>4} {2:12.6f} {3:12.6f} {4:12.6f} {5:>12.6f}".format(int(col[0]),col[1],float(col[2]),float(col[3]),float(col[4]),float(atom_type)))
                for i in col[6:(len(col))-1]:
                    outfile.write( "{0:>5} ".format(int(i)))
                outfile.write( "{0:>5}\n".format(int(col[-1])))
    print "---Done---\n"
    return

atoms = []
count_unique(template,atoms)
atoms = np.unique(atoms)
atom_no_dictionary = {}
atom_type_dictionary = {}

output_numbs = defaultdict(list)
output_symbs = defaultdict(list)
vector_dict = defaultdict(list)
file_type = template

print "Generating primary connectivities for template file...\n"

gen_1_connections(file_type,atoms,output_numbs,output_symbs,vector_dict)    
vector_dictionary(output_symbs,vector_dict)

print "Primary connectivities established\n"
if check_vector(template,vector_dict):
    print "Testing atom types...\n"
else:
    print "Atom types are unique, printing outfile\n"
    generate_dictionaries(file_type,vector_dict,atom_type_dictionary,atom_no_dictionary)
    file_type = testfile
    output_numbs = defaultdict(list)
    output_symbs = defaultdict(list)
    vector_dict = defaultdict(list)
    gen_1_connections(file_type,atoms,output_numbs,output_symbs,vector_dict)
    vector_dictionary(output_symbs,vector_dict)
    generate_dictionaries(file_type,vector_dict,atom_type_dictionary,atom_no_dictionary)
    print_output(atom_no_dictionary,atom_type_dictionary,testfile)
    exit() # Leave Code, atoms are assigned
    

connections_template = 1
while check_vector(template,vector_dict):
    connections_template += 1
    print "Atom types are not unique, generating connectivity order: {0} for template file \n...".format(connections_template)
    input_numbs = output_numbs
    input_symbs = output_symbs
    output_numbs = defaultdict(list)
    output_symbs = defaultdict(list)
    gen_n_connections(file_type,atoms,input_numbs,input_symbs,output_numbs,output_symbs,vector_dict)
    vector_dictionary(output_symbs,vector_dict)
    print "Checking uniqueness...\n"

print "Uniqueness found after {0} connectivity cycles\n".format(connections_template)
generate_dictionaries(file_type,vector_dict,atom_type_dictionary,atom_no_dictionary)
print "Generating testfile connectivities...\n"
file_type = testfile
output_numbs = defaultdict(list)
output_symbs = defaultdict(list)
vector_dict = defaultdict(list)
gen_1_connections(file_type,atoms,output_numbs,output_symbs,vector_dict)
vector_dictionary(output_symbs,vector_dict)

connections_testfile = 1
while connections_testfile < connections_template:
    connections_testfile += 1
    print "Connectivity cycle {0} for testfile\n".format(connections_testfile)
    input_numbs = output_numbs
    input_symbs = output_symbs
    output_numbs = defaultdict(list)
    output_symbs = defaultdict(list)
    gen_n_connections(file_type,atoms,input_numbs,input_symbs,output_numbs,output_symbs,vector_dict)
    vector_dictionary(output_symbs,vector_dict)


generate_dictionaries(file_type,vector_dict,atom_type_dictionary,atom_no_dictionary)
print "Printing to outfile...\n"
print_output(atom_no_dictionary,atom_type_dictionary,testfile)
    
 
       

