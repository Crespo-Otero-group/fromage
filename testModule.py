import sys
import numpy
import subprocess
import os
from atom import Atom
from readfile import *
from editinputs import *
from molselect import *

one = Atom("C",1,1,1,2)
two = Atom("C",2,1,1,-2)
three = Atom("C",1,1,2,1)

atoms = [one,two,three]
