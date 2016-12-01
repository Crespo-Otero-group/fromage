import sys
import numpy
import subprocess
import os
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

here = os.path.dirname(os.path.realpath(__file__))
popDir = "POP"
popPath = os.path.join(here,popDir)

boop = readBader(os.path.join(popPath,"ACF"),os.path.join(popPath,"naph"))
for charge in boop:
    print charge
print sum(boop)

atat = Atom("c",0,0,0,5)
print atat
a,b=atat.electrons()
print a
print b
