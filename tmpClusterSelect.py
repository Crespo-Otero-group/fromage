import sys
import numpy
import subprocess
import os
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

name = "inname"
clustRad = 10
maxBL= 1.58


vectors = readvasp(name)["vectors"]
atoms = readvasp(name)["atoms"]

N = len(atoms)
baryX, baryY, baryZ = 0, 0, 0
for atom in atoms:
    baryX += atom.x / N
    baryY += atom.y / N
    baryZ += atom.z / N
transAtoms=[]

for atom in atoms:
    transAtoms.append(atom.translate(-baryX, -baryY, -baryZ))


seedatoms=[]

for atom in transAtoms:
    if atom.dist(0, 0, 0) < clustRad:
        seedatoms.append(atom)

clustAtoms=[]

for atom in seedatoms:
    if atom not in clustAtoms:
        mol2Add = select(maxBL, transAtoms, transAtoms.index(atom))
        for atom2Add in mol2Add:
            clustAtoms.append(atom2Add)

print clustAtoms
writexyz(name+".clust",clustAtoms)
