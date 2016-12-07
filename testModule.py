import sys
import numpy
import subprocess
import os
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

name = "naphthalene222"

here = os.path.dirname(os.path.realpath(__file__))
qeDir = "QE"
qePath = os.path.join(here,qeDir)

vectors = readvasp(name)["vectors"]
atoms = readvasp(name)["atoms"]

os.chdir(qePath)


editQE("yabadabada", vectors, atoms)
