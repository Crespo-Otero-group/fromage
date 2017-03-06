import sys
import numpy
import subprocess
import os
import scipy
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

atoms=readGauss("2HC1")
for atom in atoms:
    print atom
