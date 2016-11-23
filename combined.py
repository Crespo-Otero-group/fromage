# main function which takes in a vasp unit cell file
# and outputs a Gaussian or Turbomole output file of
# the subsystem embedded in Ewald charges
import sys
import numpy
import subprocess
import os
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

##########
##########
# Start user inputs
##########
##########

# name of the job goes here
name = "naphthalene222"

# maximum bond length when defining a molecule
maxBL = 1.58

# label of an atom which will be part of the quantum cluster
# warning: [0,N-1], not [1,N], remember to take away 1 to any atom labels
labelAtom = 2

# the number of checkpoints in region 1
nChk = 1000

# the number of constrained charge atoms
# i.e. atoms in regions 1 and 2
# pick a larger number for larger "real system" atoms
nAt = 1000

# Ewald will multiply the unit cell in the direction
# of the a, b or c vector 2N times (N positive and N negative)
aN = 2
bN = 2
cN = 2

# Gaussian (0) or Turbomole (1)
program = 0

# the cluster will be of all molecules with atoms less than
# clustRad away from the centre of the central molecule
clustRad = 10

# how many times the input cluster needs to be repeated along each vector
# positively and negatively to be able to contain the cluster to select.
# the supercluster ends up being (1+2*traAN)*(1+2*traBN)*(1+2*traCN) time bigger
traAN = 1
traBN = 1
traCN = 1

##########
##########
# End user inputs
##########
##########

here = os.path.dirname(os.path.realpath(__file__))
cp2kDir = "CP2K"
cp2kPath = os.path.join(here, cp2kDir)
ewald1Dir = "EWALD1"
ewald1Path = os.path.join(here, ewald1Dir)
ewald2Dir = "EWALD2"
ewald2Path = os.path.join(here, ewald2Dir)

# extracts the cell information
# from a vasp file
vectors = readvasp(name)["vectors"]
atoms = readvasp(name)["atoms"]

# sets up a cp2k relaxation (or single point by changing the template)
editcp2k(cp2kPath, name, vectors, atoms)
# os.chdir(cp2kPath)
#subprocess.call("cp2k.popt -i cp2k."+name+".in -o cp2k."+name+".out",shell=True)
# os.chdir(here)

# reads the output of the relaxation
relaxedAtoms = readxyz(os.path.join(cp2kPath, name + "-pos-1"))[-1]
charges = readcp2k(os.path.join(cp2kPath, "cp2k." + name))["charges"]

# due to poor precision in cp2k, the molecule is usually
# electrically neutral only up to 10e-7.
# here the charge of the last atom is tweaked to compensate
if sum(charges) != 0.0:
    charges[-1] -= sum(charges)

# assigns Mulliken charges to atoms
for index, atom in enumerate(relaxedAtoms):
    atom.q = charges[index]


# the atoms which are part of the same selected molecule.
# first with intact coordinates and second with appropriate
# atoms translated to join up the molecule with their image
fullMol, fullMolTrans = select2(maxBL, relaxedAtoms, labelAtom, vectors)

# atoms from molecule which need translating
partMol = [atom for atom in fullMol if atom not in fullMolTrans]
# atoms from molecule which were translated
partMolImg = [atom for atom in fullMolTrans if atom not in fullMol]

tweakedCell = copy(relaxedAtoms)

# for all atoms in the cell
for atom in tweakedCell:
    # if the atom needs translating
    if atom in partMol:
        # translate the atom
        tweakedCell[tweakedCell.index(atom)] = partMolImg[
            partMol.index(atom)]

# now in tweakedCell if the selected molecule had been
# chopped off from the cell, the chopped off atoms
# have been translated back to form a whole molecule


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
for atom in tweakedCell:
    transAtoms.append(atom.translate(-baryX, -baryY, -baryZ))

# translate only the molecule as well
transMol = []
for atom in fullMolTrans:
    transMol.append(atom.translate(-baryX, -baryY, -baryZ))


# write Ewald input files for a small cluster

writeuc(ewald1Path, name, vectors, aN, bN, cN, transAtoms)
writeqc(ewald1Path, name, transMol)
# For now the following line ensures no defects in the
# step of the Ewald procedure
# In future versions .dc could be different to .qc
subprocess.call("cp " + os.path.join(ewald1Path, name) + ".qc " +
                os.path.join(ewald1Path, name) + ".dc", shell=True)
writeEwIn(ewald1Path, name, nChk, nAt)
writeSeed(ewald1Path)
# run Ewald
os.chdir(ewald1Path)
#subprocess.call("./Ewald < ewald.in." + name, shell=True)
os.chdir(here)
# read points output by Ewald
points = readPoints(os.path.join(ewald1Path, name))

# write Ewald input files for a large cluster

# this will contain an even bigger supercell made of 8 supercells
superMegaCell = []

traA = range(-traAN,traAN+1)
traB = range(-traBN,traBN+1)
traC = range(-traCN,traCN+1)

for i in traA:
    for j in traB:
        for k in traC:
            traX = vectors[0][0] * i + vectors[1][0] * j + vectors[2][0] * k
            traY = vectors[0][1] * i + vectors[1][1] * j + vectors[2][1] * k
            traZ = vectors[0][2] * i + vectors[1][2] * j + vectors[2][2] * k

            for atom in transAtoms:
                superMegaCell.append(atom.translate(traX, traY, traZ))




# atoms within the sphere of rad clustRad
seedatoms = []

for atom in superMegaCell:
    if atom.dist(0, 0, 0) < clustRad:
        seedatoms.append(atom)

# atoms in the cluster (seedatoms + atoms to complete molecules)
clustAtoms = []
for atom in seedatoms:
    if atom not in clustAtoms:
        mol2Add = select(maxBL, superMegaCell, superMegaCell.index(atom))
        for atom2Add in mol2Add:
            clustAtoms.append(atom2Add)


# write Ewald input files for a large cluster

writeuc(ewald2Path, name + ".clust", vectors, aN, bN, cN, transAtoms)
writeqc(ewald2Path, name + ".clust", clustAtoms)
# For now the following line ensures no defects in the
# step of the Ewald procedure
# In future versions .dc could be different to .qc
subprocess.call("cp " + os.path.join(ewald2Path, name + ".clust") + ".qc " +
                os.path.join(ewald2Path, name + ".clust") + ".dc", shell=True)
writeEwIn(ewald2Path, name + ".clust", nChk, nAt)
writeSeed(ewald2Path)
# run Ewald
os.chdir(ewald2Path)
#subprocess.call("./Ewald < ewald.in." + name + ".clust", shell=True)
os.chdir(here)
# read points output by Ewald
pointsClust = readPoints(os.path.join(ewald2Path, name + ".clust"))



outerAtoms = []

for atom in clustAtoms:
    if atom not in transAtoms:
        outerAtoms.append(atom)

outerPoints = []

for atom in outerAtoms:
    point = copy(atom)
    point.elem = "point"
    outerPoints.append(point)

pointsInner = outerPoints.append(pointsClust)

# select the program to do the high level subsystem calculation with
if program == 0:
    # write Gaussian input file
    writeGauss(name, transMol, points, 1)
    writeGauss(name + ".clust", clustAtoms, pointsClust, 0)
elif program == 1:
    # write Turbomole control
    editControl(points)
    #subprocess.call("jobex -ex -c 300",shell=True)



# at this point we have extracted the energy of the
# subsystem embedded in Ewald charges


# the next lines presume that we also want to calculate the system
# at a low level of theory in cp2k. This might turn out to be impossible
# if we want to include our own charge embedding


# sets up a cp2k calculation of the subsystem at low level theory
# currently the bounding box is 5* the input cell
# this should be plenty but it could be changed to a parameter in the future
#editcp2k(name + ".low", vectors*5, transMol)
#subprocess.call("cp2k.popt -i cp2k."+name+".low.in -o cp2k."+name+".low.out",shell=True)
#lowEnergy = readcp2k("cp2k." + name +".low")["energy"]
