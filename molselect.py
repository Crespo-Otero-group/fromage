from numpy import *
from atom import Atom
from copy import copy
# takes a maximum interatomic distance,
# list of atoms and a label as inputs
# returns a list of the atoms connected
# with the labeled atom in a molecule


def select(maxR, atoms, label):
    N = len(atoms)
    M = zeros((N, N))

    selected = [atoms[label]]
    n = True
    while n == True:
        n = False
        for i in selected:
            for j in atoms:
                if i.dist(j.x, j.y, j.z) <= maxR and j not in selected:
                    selected.append(j)
                    n = True

    return selected

# same as select but also takes in a matrix
# of lattice vectors which allows the computation
# of distances to include periodic images.
# also outputs the atoms translated
# and untranslated by lattice vectors


def select2(maxR, atoms, label, vectors):

    # list of selected atoms from the unit cell
    selected = [atoms[label]]
    # list of selected atoms where the periodic image
    # atoms are translated back to form a molecule
    selectedImg = [atoms[label]]

    # n is true as long as there are still molecules to add to the list
    n = True
    while n == True:
        n = False
        for i in selectedImg:
            for j in atoms:
                # contains the distance from the point or image and the
                # coordinates of the point or image
                gamma = i.distLat(j.x, j.y, j.z, vectors[
                                  0], vectors[1], vectors[2])

                # if the atom is close enough to be part of the molecule
                # and is not already part of the molecule
                if gamma[0] <= maxR and j not in selected:
                    selected.append(j)
                    k = copy(j)
                    k.x, k.y, k.z = gamma[1:]
                    selectedImg.append(k)
                    n = True

    return selected, selectedImg
