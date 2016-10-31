# main function which takes in a vasp unit cell file
# remember to module load cp2k
import sys
import numpy
import subprocess
from atom import Atom
from readfile import *
from editinputs import *


name="naphthalene222"

vectors = readvasp(name)["vectors"]
atoms= readvasp(name)["atoms"]

editcp2k(name, vectors,atoms)
#subprocess.call("cp2k.popt -i cp2k."+name+".in -o cp2k."+name+".out",shell=True)


relaxedAtoms = readxyz(name+"-pos-1")[-1]
charges = readcp2kMull("cp2k."+name)


for index, atom in enumerate(relaxedAtoms):
    atom.q = charges[index]
    print st
