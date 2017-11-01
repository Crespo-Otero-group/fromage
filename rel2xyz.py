#!/usr/bin/env python
import time

import sys
import argparse
import assign_charges as ac
import read_file as rf
import edit_file as ef
import handle_atoms as ha
from math import sqrt
import numpy as np

in_f = sys.argv[1]
atoms=rf.read_xyz(in_f)[-1]
vectors=np.loadtxt("vectors")
mega=ha.make_mega_cell(atoms,1,1,1,vector_mat)
ef.write_xyz(in_f[:-4]+"-mega.xyz",mega)
