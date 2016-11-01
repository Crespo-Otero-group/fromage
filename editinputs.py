# functions for editing inputs of various programs

import sys
import numpy
from atom import Atom

# edits a cp2k template file called cp2k.template.in
# and writes a new version called cp2k.[input name].in
# with required lattice vectors and atomic positions


def editcp2k(inName, vectors, atoms):
    with open("cp2k.template.in") as tempFile:
        tempContent = tempFile.readlines()

    # strings for each lattice vector
    aVec = "\t".join(map(str, vectors[0]))
    bVec = "\t".join(map(str, vectors[1]))
    cVec = "\t".join(map(str, vectors[2]))

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
                cp2kIn.write(str(atom.elem) + " \t" + str(atom.x) +
                             " \t" + str(atom.y) + " \t" + str(atom.z) + "\n")

        else:  # if no tag is found
            cp2kIn.write(line)

    tempFile.close()
    cp2kIn.close()
    return

# edits a cp2k submission script template to
# give it the correct job, input and output names
# possibly useless, stored here for later use


def editcp2kSub(inName):
    with open("script-cp2k.template") as tempFile:
        tempContent = tempFile.readlines()

    subScript = open("script-cp2k." + inName, "w")

    for line in tempContent:
        if "XXX__NAME__XXX" in line:
            subScript.write(line.replace("XXX__NAME__XXX", inName))
        else:
            subScript.write(line)

    tempFile.close()
    return

# writes an xyz file from a list of atoms


def writexyz(inName, atoms):
    outFile = open(inName + ".xyz", "w")
    outFile.write(str(len(atoms)) + "\n")
    outFile.write(inName + "\n")

    for atom in atoms:
        outFile.write(atom.xyzStr())
    return


# writes a .uc file for the Ewald program
# input the name, the matrix of lattice vectors
# each multiplication of cell through a vector
# and a list of Atom objects
def writeuc(inName, vectors, aN, bN, cN, atoms):
    outFile = open(inName + ".uc", "w")
    outFile.write("\t".join(map(str, vectors[0])) + "\t" + str(aN) + "\n")
    outFile.write("\t".join(map(str, vectors[1])) + "\t" + str(bN) + "\n")
    outFile.write("\t".join(map(str, vectors[2])) + "\t" + str(cN) + "\n")

    # Transpose to ge the transformation matrix
    M = numpy.transpose(vectors)
    # Inverse transformation matrix
    U = numpy.linalg.inv(M)

    for atom in atoms:
        dirPos = [atom.x, atom.y, atom.z]
        fracPos = numpy.dot(U, dirPos).tolist()
        for coord in fracPos:
            if coord < 0:
                fracPos[fracPos.index(coord)] = 1 + coord
        strLine = "\t".join(map(str, fracPos)) + "\t" + \
            str(atom.q) + "\t" + str(atom.elem) + "\n"
        outFile.write(strLine)

    return

# writes a .qc file for Ewald with a name and a list of atoms


def writeqc(inName, atoms):
    outFile = open(inName + ".qc", "w")
    for atom in atoms:
        outFile.write(str(atom))
