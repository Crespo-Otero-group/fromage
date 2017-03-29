# main function which takes in a vasp unit cell file
# and outputs a Gaussian output file of
# the subsystem embedded in Ewald charges
# The potential can also be calculated self-consistently
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
name = "2HC1"

# maximum bond length when defining a molecule
maxBL = 1.58

# label of an atom which will be part of the quantum cluster
# warning: [0,N-1], not [1,N]
labelAtom = 2

# the number of checkpoints in region 1
nChk = 1000

# the number of constrained charge atoms
# i.e. atoms in regions 1 and 2
nAt = 2700

# Ewald will multiply the unit cell in the direction
# of the a, b or c vector 2N times (N positive and N negative)
aN = 2
bN = 2
cN = 2

# Periodic program QE (0) CP2K (1)
programPer = 1

# Population analysis method if pertinent
# Mulliken(0) Hirshfeld(1) RESP(1)
programPop = 1

# relaxing with cp2k or just single point?
relaxBool = True

# make a vasp file at the end with an excited molecule in the unit cell?
loopBool = False

# assign excited charges and loop the Ewald calculation
ewLoop = False


# the cluster will be of all molecules with atoms less than
# clustRad away from the centre of the central molecule
clustRad = 5

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
ewaldDir = "EWALD"
ewaldPath = os.path.join(here, ewaldDir)
gaussianDir = "GAUSSIAN"
gaussianPath = os.path.join(here, gaussianDir)
popDir = "POP"
popPath = os.path.join(here, popDir)
qeDir = "QE"
qePath = os.path.join(here, qeDir)

# extracts the cell information
# from a vasp file
vectors = readvasp(name)["vectors"]
atoms = readvasp(name)["atoms"]

# sets up a periodic relaxation (or single point by changing the template and option)

# if Quantum Espresso
if programPer == 0:
    os.chdir(qePath)
    editQE(name, vectors, atoms)
    #subprocess.call("pw.x <qe."+name+".in> qe."+name+".out",shell=True)
    if relaxBool == True:
        relaxedAtoms = readQE("qe." + name)
    else:
        relaxedAtoms = atoms
    # get a population analysis
    editPP(name)
    #subprocess.call("pp.x <pp."+name+".in> pp."+name+".out",shell=True)


# if CP2K (currently don't use Bader with CP2K)
elif programPer == 1:
    os.chdir(cp2kPath)
    editcp2k(name, vectors, atoms)
    #subprocess.call("cp2k.popt -i cp2k."+name+".in -o cp2k."+name+".out",shell=True)

    # reads the output of the relaxation
    if relaxBool == True:
        relaxedAtoms = readxyz(name + "-pos-1")[-1]
    # unless they are already relaxed and it was single point
    else:
        relaxedAtoms = atoms
    charges = readcp2k("cp2k." + name,programPop)["charges"]
os.chdir(here)

# due to poor precision in cp2k, the molecule is usually
# electrically neutral only up to 10e-7.
# here the charge of the last atom is tweaked to compensate
if sum(charges) != 0.0:
    charges[-1] -= sum(charges)

# assigns charges to atoms
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


# for all atoms in the cell
for atom in relaxedAtoms:
    # if the atom needs translating
    if atom in partMol:
        # translate the atom
        relaxedAtoms[relaxedAtoms.index(atom)] = partMolImg[
            partMol.index(atom)]

# now in relaxedAtoms if the selected molecule had been
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
for atom in relaxedAtoms:
    transAtoms.append(atom.translate(-baryX, -baryY, -baryZ))

# translate only the molecule as well
transMol = []
for atom in fullMolTrans:
    transMol.append(atom.translate(-baryX, -baryY, -baryZ))





# this will contain an even bigger supercell
superMegaCell = []

traA = range(-traAN,traAN+1)
traB = range(-traBN,traBN+1)
traC = range(-traCN,traCN+1)

# multiplying the cell
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


clustPoints = []
shellAtoms= []
for atom in clustAtoms:
    if atom not in transMol:
        shellAtoms.append(atom)
        newClustPoint=copy(atom)
        newClustPoint.elem=("point")
        clustPoints.append(newClustPoint)

os.chdir(here)
writexyz("mol", transMol)
writexyz("clust", clustAtoms)
writexyz("shell", shellAtoms)


###EWALD

os.chdir(ewaldPath)
writeuc(name, vectors, aN, bN, cN, transAtoms)
writeqc(name, clustAtoms)
# For now the following line ensures no defects in the
# step of the Ewald procedure
# In future versions .dc could be different to .qc
subprocess.call("cp " + name + ".qc " + name + ".dc", shell=True)
writeEwIn(name, nChk, nAt)
writeSeed()
# run Ewald
subprocess.call("./Ewald < ewald.in." + name, shell=True)
# read points output by Ewald
points = readPoints(name)
os.chdir(here)


###GAUSSIAN


# write Gaussian input file
os.chdir(gaussianPath)
writeGauss(name+".clust", clustAtoms, points)
writeGauss(name, transMol, points+clustPoints)
#subprocess.call("g09 "+name+".com",shell=True)
#subprocess.call("formchk "+name+".chk "+name+".fck",shell=True)
#subprocess.call("cubegen 1 Density=CI "+name+".fck "+name+".cube",shell=True)
