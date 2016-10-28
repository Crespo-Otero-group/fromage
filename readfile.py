# functions to read files of different formats

import sys
import numpy
from atom import Atom


# Takes the prefix of a vasp file and returns a matrix of lattice vectors
# and atoms as a list of Atom objects
def readvasp(inName):
    atoms = []
    with open(inName + ".vasp") as vaspFile:
        vaspContent = vaspFile.readlines()

    # Lattice vectors
    # Make sure the vasp "lattice constant" scaling is set to 1.0
    vec1 = vaspContent[2].split()
    vec2 = vaspContent[3].split()
    vec3 = vaspContent[4].split()

    # Matrix from vectors
    M = numpy.zeros((3, 3))
    M[0] = vec1
    M[1] = vec2
    M[2] = vec3

    # Reads names of elements and amounts
    species = vaspContent[5].split()
    amountsStr = vaspContent[6].split()
    amounts = map(int, amountsStr)

    # Make Atom objects from file

    for element in species:

        # Position of the first and last atom of one type
        # in the vasp file
        firstAt = 8 + sum(amounts[:species.index(element)])
        lastAt = 8 + sum(amounts[:species.index(element) + 1])

        for line in vaspContent:
            if vaspContent.index(line) in range(firstAt, lastAt):
                xAtom, yAtom, zAtom = line.split()
                atoms.append(Atom(element, xAtom, yAtom, zAtom))

    return {"vectors": M, "atoms": atoms}
