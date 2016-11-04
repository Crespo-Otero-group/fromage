# functions for editing and writing inputs of various programs

import sys
import numpy
from atom import Atom
from random import randint

# edits a cp2k template file called cp2k.template.in
# and writes a new version called cp2k.[input name].in
# with required lattice vectors and atomic positions


def editcp2k(inName, vectors, atoms):
    with open("cp2k." + inName + ".template.in") as tempFile:
        tempContent = tempFile.readlines()

    # strings for each lattice vector
    aVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[0])
    bVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[1])
    cVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[2])

    cp2kIn = open("cp2k." + inName + ".in", "w")

    for line in tempContent:

        # writes the name of the calculation at the top of the file
        if "XXX__NAME__XXX" in line:
            cp2kIn.write(line.replace("XXX__NAME__XXX", inName))

        # replace the tags with the coordinates of lattice vectors
        elif "XXX__AVEC__XXX" in line:
            cp2kIn.write(line.replace("XXX__AVEC__XXX", aVec))
        elif "XXX__BVEC__XXX" in line:
            cp2kIn.write(line.replace("XXX__BVEC__XXX", bVec))
        elif "XXX__CVEC__XXX" in line:
            cp2kIn.write(line.replace("XXX__CVEC__XXX", cVec))

        # writes atomic coordinates
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                lineStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z)
                cp2kIn.write(lineStr + "\n")

        else:  # if no tag is found
            cp2kIn.write(line)

    cp2kIn.close()
    return


# writes an xyz file from a list of atoms
# for future use


def writexyz(inName, atoms):
    outFile = open(inName + ".xyz", "w")
    outFile.write(str(len(atoms)) + "\n")
    outFile.write(inName + "\n")

    for atom in atoms:
        outFile.write(atom.xyzStr() + "\n")
    outFile.close()
    return


# writes a .uc file for the Ewald program
# input the name, the matrix of lattice vectors
# each multiplication of cell through a vector
# and a list of Atom objects


def writeuc(inName, vectors, aN, bN, cN, atoms):

    line1 = vectors[0].tolist() + [aN]
    line2 = vectors[1].tolist() + [bN]
    line3 = vectors[2].tolist() + [cN]
    outFile = open(inName + ".uc", "w")
    outFile.write("{:10.6f} {:10.6f} {:10.6f} {:10d}".format(*line1) + "\n")
    outFile.write("{:10.6f} {:10.6f} {:10.6f} {:10d}".format(*line2) + "\n")
    outFile.write("{:10.6f} {:10.6f} {:10.6f} {:10d}".format(*line3) + "\n")

    # transpose to get the transformation matrix
    M = numpy.transpose(vectors)
    # inverse transformation matrix
    U = numpy.linalg.inv(M)

    for atom in atoms:
        # change of basis transformation
        dirPos = [atom.x, atom.y, atom.z]
        fracPos = numpy.dot(U, dirPos).tolist()
        for coord in fracPos:
            # if the coordinate is negative
            if coord < 0:
                # translate it to the range [0,1]
                fracPos[fracPos.index(coord)] = 1 + coord
        strLine = "{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:>6}".format(
            *fracPos + [atom.q] + [atom.elem]) + "\n"
        outFile.write(strLine)
    outFile.close()
    return


# writes a .qc file for Ewald with a name and a list of atoms


def writeqc(inName, atoms):
    outFile = open(inName + ".qc", "w")
    for atom in atoms:
        outFile.write(str(atom) + "\n")
    outFile.close()
    return


# writes a ewald.in file from the job name,
# the amount of checkpoints in zone 1 and
# the amount of atoms with constrained charge


def writeEwIn(inName, nChk, nAt):
    outFile = open("ewald.in." + inName, "w")
    outFile.write(inName + "\n")
    outFile.write(str(nChk) + "\n")
    outFile.write(str(nAt) + "\n")
    outFile.write("0\n")
    outFile.close()
    return

# writes a seed file for Ewald


def writeSeed():
    outFile = open("seedfile", "w")
    seed1 = randint(1, 2**31 - 86)
    seed2 = randint(1, 2**31 - 250)
    outFile.write(str(seed1) + " " + str(seed2))
    outFile.close()


# writes a Gaussian input file from a template
# with atoms and point charges as inputs


def writeGauss(inName, atoms, points):
    with open("template.com") as tempFile:
        tempContent = tempFile.readlines()

    outFile = open(inName + ".com", "w")

    for line in tempContent:
        if "XXX__NAME__XXX" in line:
            outFile.write(line.replace("XXX__NAME__XXX", inName))
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                # atomStr = str(atom.elem) + " \t" + str(atom.x) + \
                #    " \t" + str(atom.y) + " \t" + str(atom.z) + "\n"
                outFile.write(atomStr)
        elif "XXX__CHARGES__XXX" in line:
            for point in points:
                pointStr = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    point.x, point.y, point.z, point.q) + "\n"
                # pointStr = str(point.x) + " \t" + str(point.y) + \
                #    " \t" + str(point.z) + " \t" + str(point.q) + "\n"
                outFile.write(pointStr)
        else:
            outFile.write(line)
    outFile.close()
    return


def editControl(points):
    with open("control.template") as ctrlFile:
        content = ctrlFile.readlines()

    outFile = open("control", "w")

    for line in content:
        if "$end" in line:
            outFile.write("$point_charges\n")
            for point in points:
                pointStr = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    point.x, point.y, point.z, point.q) + "\n"
                outFile.write(pointStr)
            outFile.write("$end")
        else:
            outFile.write(line)

    return
