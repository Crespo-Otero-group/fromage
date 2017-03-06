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

    # lattice vectors

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
                xAtom, yAtom, zAtom = map(float, line.split())
                atoms.append(Atom(element, xAtom, yAtom, zAtom))
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
        if line:
            if line.split()[0].isdigit():

                # list of atom objects inside on relaxation step
                atoms = []

                # from 2 lines after the amount of atoms to the last atom line
                # for the relaxation step
                for lineInStep in xyzContent[xyzContent.index(line) + 2:xyzContent.index(line) + int(line) + 2]:
                    elemAtom = lineInStep.split()[0]
                    xAtom, yAtom, zAtom = map(float, lineInStep.split()[1:])
                    atoms.append(Atom(elemAtom, xAtom, yAtom, zAtom))

                atomStep.append(atoms)

    xyzFile.close()
    return atomStep


# reads a cp2k output file to extract the Mulliken charges in a list


def readcp2k(inName):
    with open(inName + ".out") as cp2kFile:
        cp2kContent = cp2kFile.readlines()

    # find last occurrence of Mulliken charges
    lastMull = len(cp2kContent) - 1 - \
        cp2kContent[::-1].index("                           Hirshfeld Charges\n")

    charges = []

    for line in cp2kContent[lastMull + 3:]:
        if line != '\n':
            charges.append(float(line.split()[5]))
        else:
            break
    # find each occurrence of Energy
    for line in cp2kContent:
        if "ENERGY|" in line:
            energy = float(line.split()[8])

    cp2kFile.close()
    return {"charges": charges, "energy": energy}


# returns point charges as Atom objects of element
# "point". the format is .pts-tb, output from the house
# version of Ewald


def readPoints(inName):
    with open(inName + ".pts-tb") as ptsFile:
        ptsContent = ptsFile.readlines()

        # store point charges here
    points = []

    for line in ptsContent:
        xIn, yIn, zIn, qIn = map(float, line.split())
        point = Atom("point", xIn, yIn, zIn, qIn)
        points.append(point)

    return points


# Same as readcp2k but for gaussian log files

def readGMull(inName):
    with open(inName + ".log") as gaussFile:
        content = gaussFile.readlines()

    # find last occurrence of Mulliken charges
    lastMull = len(content) - 1 - \
        content[::-1].index(" Mulliken charges:\n")

    charges=[]

    for line in content[lastMull+2:]:
        if line.split()[0].isdigit():
            charges.append(float(line.split()[2]))
        else:
            break
    # find each occurrence of Energy
    for line in content:
        if "Total Energy" in line:
            energy = float(line.split()[4])

    return {"charges": charges, "energy": energy}

# returns the charge values from a bader calculation

def readBader(inDat):
    with open(inDat + ".dat") as baderFile:
        contentB = baderFile.readlines()

    # electron charge per atom
    charges = []
    for line in contentB:
        if line.split()[0].isdigit():
            charge = float(line.split()[4])
            charges.append(charge)

    return charges


# returns final positions of a Quantum Espresso calculation
def readQE(inFile):
    with open(inFile+".out") as fileQE:
        content = fileQE.readlines()

    lastPos = 0
    for line in content[::-1]:
        if "ATOMIC_POSITIONS" in line.split():
            lastPos = content[::-1].index(line)
            break

    atoms = []
    for line in content[-lastPos:]:
        if line == "End final coordinates\n":
            break
        elem, xPos, yPos, zPos = line.split()
        atom2Add = Atom(elem, xPos, yPos, zPos, 0)
        atoms.append(atom2Add)
    return atoms

# returns the atom list in a Gaussian input file
def readGauss(inFile):
    with open(inFile+".com") as fileGauss:
        content= fileGauss.readlines()

    for line in content:
        if line != "\n":
            if line.split()[0].isdigit():
                lastPos=content.index(line)
                break
    atoms = []
    for line in content[lastPos+1:]:
        if line =="\n" or not line:
            break
        elem, xPos, yPos, zPos = line.split()
        atom2Add = Atom(elem, xPos, yPos, zPos, 0)
        atoms.append(atom2Add)
    return atoms
