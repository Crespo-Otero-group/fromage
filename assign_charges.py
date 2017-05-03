from atom import Atom
from readfile import *
from collections import Counter
import numpy as np
from copy import deepcopy

# takes in a list of atoms, some lattice vectors and a max bond length
# returns a matrix of first connectivities
def detect_1_connect(in_atoms, vectors, max_BL):
    nat_mol = len(in_atoms)
    cnct = np.zeros((nat_mol, nat_mol))
    for i, i_atom in enumerate(in_atoms):
        for j, j_atom in enumerate(in_atoms):
            if vectors is None:
                if i_atom.dist(j_atom.x, j_atom.y, j_atom.z) < max_BL:
                    cnct[i][j] = 1
            else:
                if i_atom.distLat(j_atom.x, j_atom.y, j_atom.z, vectors[0], vectors[1], vectors[2])[0] < max_BL:
                    cnct[i][j] = 1
    return cnct

# returns connectivity matrix of one higher order than before


def expand_connect(in_mat):
    out_mat = np.copy(in_mat)
    for i, row in enumerate(in_mat):
        # indeces of unconnected atoms
        dangles = []
        # indeces of connected atoms
        connectors = []
        for j, element in enumerate(row):
            if element == 0 and j > i:
                dangles.append(j)
            elif element != 0:
                connectors.append(j)
        for dangle in dangles:
            orders = []
            for k, dangle_element in enumerate(in_mat[dangle]):
                if dangle_element != 0 and k in connectors:
                    orders.append(dangle_element + in_mat[i][k])
                # if orders is not empty
                if orders:
                    out_mat[i][dangle] = min(orders)
                    out_mat[dangle][i] = min(orders)
    return out_mat

# expands the connectivity matrix until it is complete
# returns the complete matrix


def complete_expand(in_mat):
    mat = np.copy(in_mat)
    i = 1
    while True:
        i += 1
        temp_mat = expand_connect(mat)
        if np.array_equal(mat, temp_mat):
            break
        mat = np.copy(temp_mat)

    return mat

# reads a list of atoms and of unique kinds and assigns the average charge
# of the atoms of a specific kind to a tuple (charge,kind)

def charged_kinds(in_atoms,in_kinds):
    q_kinds=[]
    for kind in in_kinds:
        charges=[]
        for atom in in_atoms:
            if atom.kind==kind:
                charges.append(atom.q)
        if charges: #if not empty
            avg_charge=sum(charges)/float(len(charges))
        else:
            avg_charge=0
        q_kinds.append((avg_charge,kind))
    return q_kinds

# assigns the charges of one first molecule to a set of uncharged molecules
# based on complete connectivity
# requires a maximum bond length for the connectivity definition
# use vectors = None if the atoms are not in a periodic environment
def assign_charges(char_atoms,char_vectors,unchar_atoms,unchar_vectors,bl):
    #detect the charged atom's connectivity matrix
    char_first = detect_1_connect(char_atoms, char_vectors, bl)
    char_cnct = complete_expand(char_first)

    #get charged atom kinds as a result
    kinds=[]
    for i,atom in enumerate(char_atoms):
        atom.set_connectivity(char_atoms, char_cnct[i])
        kinds.append(atom.kind)
    kinds = set(kinds)
    q_kinds=charged_kinds(char_atoms,kinds)


    #detect uncharged atom connectivity
    unchar_first=detect_1_connect(unchar_atoms,unchar_vectors,bl)
    unchar_cnct=complete_expand(unchar_first)

    #determine kind and cross check with charged kinds
    for i,atom in enumerate(unchar_atoms):
        atom.set_connectivity(unchar_atoms, unchar_cnct[i])
        for q_kind in q_kinds:
            if atom.kind==q_kind[1]:
                atom.q=q_kind[0]

    return unchar_atoms
# mol_atoms = readxyz("mol")[0]
#
# charges = readGMull("H_hf321")["charges"]
# for i, atom in enumerate(mol_atoms):
#     atom.q = charges[i]
# vectors = readvasp("2HC1")["vectors"]
# naked_atoms = readvasp("2HC1")["atoms"]
#
#
# bl = 1.9
# # numpy.set_printoptions(threshold='nan')
# mol_first = detect_1_connect(mol_atoms, None, bl)
# mol_cnct = complete_expand(mol_first)
#
# #get single mol kinds
#
# kinds=[]
# for i,atom in enumerate(mol_atoms):
#     atom.set_connectivity(mol_atoms, mol_cnct[i])
#     kinds.append(atom.kind)
# kinds = set(kinds)
#
# q_kinds=charged_kinds(mol_atoms,kinds)
#
#
# naked_first=detect_1_connect(naked_atoms,vectors,bl)
# naked_cnct=complete_expand(naked_first)
#
# # get uc kinds
#
# for i,atom in enumerate(naked_atoms):
#     atom.set_connectivity(naked_atoms, naked_cnct[i])
#     for q_kind in q_kinds:
#         if atom.kind==q_kind[1]:
#             atom.q=q_kind[0]
#
# print naked_atoms
# print mol_atoms
