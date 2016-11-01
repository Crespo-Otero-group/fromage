from numpy import *
from atom import Atom

# takes a maximum interatomic distance,
# list of atoms and a label as inputs
# returns a list of the atoms connected
# with the labeled atom in a molecule


def select(maxR, atoms, label):
    N = len(atoms)
    M = zeros((N, N))

    selected=[atoms[label]]
    n=0
    while n==0:
        n=1
        for i in selected:
            for j in atoms:
                if i.dist(j.x, j.y, j.z) <= maxR and j not in selected:
                    selected.append(j)
                    n=0

    print selected
    print len(selected)
    return

# same as select but also takes in a matrix
# of lattice vectors which allows the computation
# of distances to include periodic images

def select2(maxR, atoms, label, vectors):
    N = len(atoms)
    M = zeros((N, N))

    selected=[atoms[label]]
    n=0
    while n==0:
        n=1
        for i in selected:
            for j in atoms:
                if i.distLat(j.x, j.y, j.z, vectors[0], vectors[1], vectors[2]) <= maxR and j not in selected:
                    selected.append(j)
                    n=0

    print selected
    print len(selected)
    return
