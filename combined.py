#main function which takes in a vasp unit cell file

import sys
import numpy
from atom import Atom
from readfile import *
from editinputs import *


name="naphthalene222"

vectors = readvasp(name)["vectors"]
atoms= readvasp(name)["atoms"]
print "hello"

editcp2k(name, vectors,atoms)
