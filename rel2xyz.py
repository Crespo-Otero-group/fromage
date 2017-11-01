#!/usr/bin/env python
import time

import sys
import read_file as rf
import edit_file as ef
import numpy as np

QE_file = sys.argv[1]
atoms=rf.read_qe(QE_file)
ef.write_xyz(in_f[:-4]+"-opt.xyz",atoms)

vectors=np.loadtxt("vectors")
mega=ha.make_mega_cell(atoms,1,1,1,vector_mat)
ef.write_xyz(in_f[:-4]+"-mega.xyz",mega)
