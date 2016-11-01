# main function which takes in a vasp unit cell file
# remember to module load cp2k
import sys
import numpy
import subprocess
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

# the name of the job goes here
name = "naphthalene222"

# extracts the unrelaxed cell information
# from a vasp file
vectors = readvasp(name)["vectors"]
atoms = readvasp(name)["atoms"]

# sets up a cp2k relaxation
editcp2k(name, vectors, atoms)
#subprocess.call("cp2k.popt -i cp2k."+name+".in -o cp2k."+name+".out",shell=True)

# reads the output of the relaxation
relaxedAtoms = readxyz(name + "-pos-1")[-1]
charges = readcp2kMull("cp2k." + name)


for index, atom in enumerate(relaxedAtoms):
    atom.q = charges[index]


# the atoms which are part of the same selected molecule
# first with intact coordiantes and second with certain
# atoms translated to join up the molecule with their image
fullMol, fullMolTrans = select2(2.2, relaxedAtoms, 2, vectors)

# atoms from molecule which need translating
partMol = [atom for atom in fullMol if atom not in fullMolTrans]
# atoms from molecule which were translated
partMolImg = [atom for atom in fullMolTrans if atom not in fullMol]

print fullMol
print fullMolTrans
print partMol
print partMolImg

# for all atoms in the cell
for atom in relaxedAtoms:
    # if the atom needs translating
    if atom in partMol:
        # translate the atom
        relaxedAtoms[relaxedAtoms.index(atom)] = partMolImg[
            partMol.index(atom)]

# finding the barycentre of the complete molecule
N = len(fullMolTrans)
baryX, baryY, baryZ = 0, 0, 0
for atom in fullMolTrans:
    baryX += atom.x / N
    baryY += atom.y / N
    baryZ += atom.z / N

# translate the whole system to have the barycentre
# at the origin
transAtoms = []
for atom in relaxedAtoms:
    transAtoms.append(atom.translate(-baryX, -baryY, -baryZ))

# translate only the molecule as well
transMol = []
for atom in fullMolTrans:
    transMol.append(atom.translate(-baryX, -baryY, -baryZ))

# write Ewald input files

# Ewald will multiply the unit cell in the direction
# of the a, b or c vector 2N times (positive and negative)
aN = 2
bN = 2
cN = 2

writeuc(name, vectors, aN, bN, cN, transAtoms)
writeqc(name, transMol)
