#main function which takes in a vasp unit cell file

import sys
import numpy
from atom import Atom
from readfile import *

print readvasp("naphthalene222")["vectors"]
