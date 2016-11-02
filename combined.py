# main function which takes in a vasp unit cell file
# remember to module load cp2k
import sys
import numpy
import subprocess
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

##########
##########
# Start user inputs
##########
##########

# the name of the job goes here
name = "naphthalene222"

# maximum bond length when defining a molecule
maxBL = 2.2

# label of an atom which will be part of the quantum cluster
labelAtom = 2

# the number of checkpoints in region 1
nChk = 1000

# the number of constrained charge atoms
# i.e. atoms in regions 1 and 2
nAt = 500

# Ewald will multiply the unit cell in the direction
# of the a, b or c vector 2N times (N positive and N negative)
aN = 2
bN = 2
cN = 2

##########
##########
# End user inputs
##########
##########



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

# due to poor precision in cp2k, the molecule is usually
# electrically neutral only up to 10e-7.
# here the charge of the last atom is tweaked to compensate
if sum(charges) != 0.0:
    charges[-1] -= sum(charges)

for index, atom in enumerate(relaxedAtoms):
    atom.q = charges[index]


# the atoms which are part of the same selected molecule
# first with intact coordiantes and second with certain
# atoms translated to join up the molecule with their image
fullMol, fullMolTrans = select2(maxBL, relaxedAtoms, labelAtom, vectors)

# atoms from molecule which need translating
partMol = [atom for atom in fullMol if atom not in fullMolTrans]
# atoms from molecule which were translated
partMolImg = [atom for atom in fullMolTrans if atom not in fullMol]



# for all atoms in the cell
for atom in relaxedAtoms:
    # if the atom needs translating
    if atom in partMol:
        # translate the atom
        relaxedAtoms[relaxedAtoms.index(atom)] = partMolImg[
            partMol.index(atom)]

# now in relaxedAtoms if the selected molecule had been
# chopped off from the cell, the chopped part is restored
# to the site of the molecule

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

writexyz("testnaph", transAtoms)
# write Ewald input files

writeuc(name, vectors, aN, bN, cN, transAtoms)
writeqc(name, transMol)
writeEwIn(name, nChk, nAt)
writeSeed()
# run Ewald
#subprocess.call("./Ewald < ewald.in."+name,shell=True)

# read points output by Ewald
points = readPoints(name)

# write Gaussian input file
writeGauss(name, transMol, points)
