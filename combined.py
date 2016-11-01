# main function which takes in a vasp unit cell file
# remember to module load cp2k
import sys
import numpy
import subprocess
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *


name="naphthalene222"
vectors = readvasp(name)["vectors"]
atoms= readvasp(name)["atoms"]


editcp2k(name, vectors,atoms)
#subprocess.call("cp2k.popt -i cp2k."+name+".in -o cp2k."+name+".out",shell=True)



relaxedAtoms = readxyz(name+"-pos-1")[-1]
charges = readcp2kMull("cp2k."+name)


for index, atom in enumerate(relaxedAtoms):
    atom.q = charges[index]


fullMol, fullMolTrans = select2(2.2, relaxedAtoms, 2, vectors)
partMol = [atom for atom in fullMol if atom not in fullMolTrans]
partMolImg = [atom for atom in fullMolTrans if atom not in fullMol]

print relaxedAtoms

for atom in relaxedAtoms:
    if atom in partMol:
        relaxedAtoms[relaxedAtoms.index(atom)] = partMolImg[partMol.index(atom)]

N=len(fullMolTrans)
baryX, baryY, baryZ = 0, 0, 0
for atom in fullMolTrans:
    baryX += atom.x/N
    baryY += atom.y/N
    baryZ += atom.z/N
print baryX, baryY, baryZ

transAtoms=[]
for atom in relaxedAtoms:
    transAtoms.append(atom.translate(-baryX, -baryY, -baryZ))

transMol=[]
for atom in fullMolTrans:
    transMol.append(atom.translate(-baryX, -baryY, -baryZ))
writexyz("testxyz", transMol)
