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
maxBL = 2.2

# label of an atom which will be part of the quantum cluster
# warning: [0,N-1], not [1,N]
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

# Periodic program QE (0) CP2K (1)
programPer = 0

# Excited state program Gaussian (0) or Turbomole (1)
programEx = 0

# Population analysis method if pertinent
# none(0) Bader(1) Hirshfeld (2)
programPop = 1

# relaxing with cp2k or just single point?
relaxBool = False

# make a vasp file at the end with an excited molecule in the unit cell?
loopBool = False

##########
##########
# End user inputs
##########
##########

here = os.path.dirname(os.path.realpath(__file__))
cp2kDir = "CP2K"
cp2kPath = os.path.join(here,cp2kDir)
ewaldDir = "EWALD"
ewaldPath = os.path.join(here,ewaldDir)
gaussianDir = "GAUSSIAN"
gaussianPath = os.path.join(here,gaussianDir)
popDir = "POP"
popPath = os.path.join(here,popDir)
qeDir = "QE"
qePath = os.path.join(here,qeDir)

# extracts the cell information
# from a vasp file
vectors = readvasp(name)["vectors"]
atoms = readvasp(name)["atoms"]

# sets up a periodic relaxation (or single point by changing the template)

# if Quantum Espresso
if programPer == 0:
    os.chdir(qePath)
    editQE(name, vectors, atoms)
    #subprocess.call("pw.x <qe."+name+".in> qe."+name+".out",shell=True)
    if relaxBool ==True:
        relaxedAtoms = readQE("qe."+name)
    else:
        relaxedAtoms = atoms
    # get a population analysis
    editPP(name)
    #subprocess.call("pp.x <pp."+name+".in> pp."+name+".out",shell=True)




# if CP2K
elif programPer == 1:
    os.chdir(cp2kPath)
    editcp2k(name, vectors, atoms)
    #subprocess.call("cp2k.popt -i cp2k."+name+".in -o cp2k."+name+".out",shell=True)


    # reads the output of the relaxation
    if relaxBool ==True:
        relaxedAtoms = readxyz(name + "-pos-1")[-1]
    # unless they are already relaxed and it was single point
    else:
        relaxedAtoms = atoms
    charges = readcp2k("cp2k." + name)["charges"]
os.chdir(here)

#if Bader
if programPop == 1:
    os.chdir(popPath)
    #subprocess.call("bader -vac off "+os.path.join(qePath,name+".cube"),shell=True)
    chargesV = readBader("ACF")
    charges = [round(b.electrons()[0] - a,6) for a, b in zip(chargesV,relaxedAtoms)]
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


# write Ewald input files

os.chdir(ewaldPath)
writeuc(name, vectors, aN, bN, cN, transAtoms)
writeqc(name, transMol)
# For now the following line ensures no defects in the
# step of the Ewald procedure
# In future versions .dc could be different to .qc
subprocess.call("cp " + name + ".qc " + name + ".dc", shell=True)
writeEwIn(name, nChk, nAt)
writeSeed()
# run Ewald
subprocess.call("./Ewald < ewald.in." + name, shell=True)
# read points output by Ewald
points = readPoints(os.path.join(name))
os.chdir(here)



# select the program to do the high level subsystem calculation with
if programEx == 0:
    # write Gaussian input file
    os.chdir(gaussianPath)
    writeGauss(name, transMol, points)
    #subprocess.call("g09 "+name+".com",shell=True)
    os.chdir(here)

elif programEx == 1:
    # write Turbomole control
    editControl(points)
    #subprocess.call("jobex -ex -c 300",shell=True)


# if we want the self-consistent version by Wilbraham
if loopBool==True:
    os.chdir(gaussianPath)
    # convert Gaussian output to xyz
    subprocess.call("g092xyz.pl "+ name+".log",shell=True)
    excitedAtoms = readxyz("g09-result")[-1]
    os.chdir(here)

    # this will be the translated cell but the molecule at the origin will be
    # relaxed in excited state
    newCell = []
    for atom in transAtoms:
        # if the atom from the cell was not part of the main molecule
        if atom not in transMol:
            # add the atom
            newCell.append(atom)
        # if the atom is in the main molecule
        else:
            # add the excited version of it instead
            atom2append=excitedAtoms[transMol.index(atom)]
            newCell.append(atom2append)
    # make a new vasp file with the information
    editVaspPos(name, newCell)
