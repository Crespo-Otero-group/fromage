# functions to read files of different formats

import sys
import numpy
from atom import Atom


# takes the prefix of a vasp file and returns a matrix of lattice vectors
# and atoms as a list of Atom objects
def readvasp(inName):

    with open(inName + ".vasp") as vaspFile:
        vaspContent = vaspFile.readlines()

    # make sure the vasp "lattice constant" scaling is set to 1.0
    # selective dynamics is not enabled
    # and the file is in Cartesian coordinates

    # sattice vectors

    vec1 = vaspContent[2].split()
    vec2 = vaspContent[3].split()
    vec3 = vaspContent[4].split()

    # matrix from vectors
    M = numpy.zeros((3, 3))
    M[0] = vec1
    M[1] = vec2
    M[2] = vec3

    # reads names of elements and amounts
    species = vaspContent[5].split()
    amountsStr = vaspContent[6].split()
    amounts = map(int, amountsStr)

    # make Atom objects from file
    atoms = []
    for element in species:

        # position of the first and last atom of one type
        # in the vasp file
        firstAt = 8 + sum(amounts[:species.index(element)])
        lastAt = 8 + sum(amounts[:species.index(element) + 1])

        for line in vaspContent:
            if vaspContent.index(line) in range(firstAt, lastAt):
                xAtom, yAtom, zAtom = line.split()
                atoms.append(Atom(element, xAtom, yAtom, zAtom))
    vaspFile.close()
    return {"vectors": M, "atoms": atoms}


# takes the prefix of an xyz file and returns
# an array of subarrayz of atoms where each
# subarray represents the system at a relaxation
# step
def readxyz(inName):

    with open(inName + ".xyz") as xyzFile:
        xyzContent = xyzFile.readlines()

    # main list where each element is a relaxation step
    atomStep = []

    for line in xyzContent:

        # if the line is the amount of atoms in the system
        if line.split()[0].isdigit():

            # list of atom objects inside on relaxation step
            atoms = []

            # from 2 lines after the amount of atoms to the last atom line
            # for the relaxation step
            for lineInStep in xyzContent[xyzContent.index(line) + 2:xyzContent.index(line) + int(line) + 2]:
                elemAtom, xAtom, yAtom, zAtom = lineInStep.split()
                atoms.append(Atom(elemAtom, xAtom, yAtom, zAtom))

            atomStep.append(atoms)

    xyzFile.close()
    return atomStep

# reads a cp2k output file to extract the Mulliken charges in a list


def readcp2kMull(inName):

    with open(inName + ".out") as cp2kFile:
        cp2kContent = cp2kFile.readlines()

    # find last occurrence of Mulliken charges
    lastMull = len(cp2kContent) - 1 - \
        cp2kContent[::-1].index(" MULLIKEN POPULATION ANALYSIS\n")

    charges = []

    for line in cp2kContent[lastMull + 3:]:
        if line.split()[0] != '#':
            charges.append(line.split()[4])
        else:
            break

    cp2kFile.close()
    return charges
