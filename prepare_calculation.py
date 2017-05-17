"""This module is a tool to prepare template and input files for cryspy.

You will need:
- config file
- .xyz file of a unit cell
- cp2k output file
- .template files for ml, mh, mg, rl calculations

And receive
- .temp files for ml, mh, mg, rl
"""
import subprocess
import os
import numpy as np
import read_file as rf
import edit_file as ef
import handle_atoms as ha
from atom import Atom
from datetime import datetime

# print start time
start_time = datetime.now()
print("STARTING TIME: " + str(start_time))

# directories
here = os.path.dirname(os.path.realpath(__file__))
ewad_dir = "EWALD"
ewald_path = os.path.join(here, ewad_dir)

# read config inputs
inputs = rf.read_config("config")

# name of the job goes here
name = inputs["name"]

# lattice vectors
a_vec = inputs["a_vec"]
b_vec = inputs["b_vec"]
c_vec = inputs["c_vec"]

vectors = np.zeros((3, 3))
vectors[0] = a_vec
vectors[1] = b_vec
vectors[2] = c_vec

print "Vectors read in config:"
print vectors
# name of the cell xyz file
if "cell_file" in inputs:
    cell_file = inputs["cell_file"]
else:
    cell_file = name + "xyz"

# name of the cp2k file with population information
if "cp2k_file" in inputs:
    cp2k_file = inputs["cp2k_file"]
else:
    cp2k_file = name + "out"

# maximum bond length when defining a molecule
if "max_bl" in inputs:
    max_bl = float(inputs["max_bl"])
else:
    max_bl = 1.7

# label of an atom which will be part of the quantum cluster
# warning: [0,N-1], not [1,N]
if "label_atom" in inputs:
    label_atom = int(inputs["label_atom"])
else:
    label_atom = 0

# the number of checkpoints in region 1
if "nchk" in inputs:
    nChk = int(inputs["nchk"])
else:
    nChk = 1000

# the number of constrained charge atoms
# i.e. atoms in regions 1 and 2
if "nat" in inputs:
    nAt = int(inputs["nat"])
else:
    nAt = 500

# Ewald will multiply the unit cell in the direction
# of the a, b or c vector 2N times (N positive and N negative)
if "an" in inputs:
    aN = int(inputs["an"])
else:
    aN = 2
if "bn" in inputs:
    bN = int(inputs["bn"])
else:
    bN = 2
if "cn" in inputs:
    cN = int(inputs["cn"])
else:
    cN = 2

# Population analysis method if pertinent
# Mulliken(0) Hirshfeld(1) RESP(2)
if "pop_method" in inputs:
    pop_method = int(inputs["pop_method"])
else:
    pop_method = 0

# the cluster will be of all molecules with atoms less than
# clust_rad away from the centre of the central molecule
if "clust_rad" in inputs:
    clust_rad = float(inputs["clust_rad"])
else:
    clust_rad = 5

# how many times the input cluster needs to be repeated along each vector
# positively and negatively to be able to contain the cluster to select.
# the supercluster ends up being (1+2*traAN)*(1+2*traBN)*(1+2*traCN) times
# bigger
if "traan" in inputs:
    traAN = int(inputs["traan"])
else:
    traAN = 2
if "trabn" in inputs:
    traBN = int(inputs["trabn"])
else:
    traBN = 2
if "tracn" in inputs:
    traCN = int(inputs["tracn"])
else:
    traCN = 2


# end config inputs

# read the input atoms
atoms = rf.read_xyz(cell_file)[-1]
print("Read " + str(len(atoms)) + " atoms in cell_file")
# read charges
charges = rf.read_cp2k(cp2k_file, pop_method)[0]
print("Read " + str(len(atoms)) + " charges in cp2k_file")

# in case there are more charges than atoms for some reason
charges = charges[:len(atoms)]
# correct charges if they are not perfectly neutral
if sum(charges) != 0.0:
    print("Charge correction: " + str(sum(charges)))
    charges[-1] -= sum(charges)

# assigns charges to atoms
for index, atom in enumerate(atoms):
    atom.q = charges[index]

# the molecule of interest and the atoms which now contain
# the full, unchopped molecule
# NB: all objects in mol are also referenced inside atoms
mol, atoms = ha.complete_mol(max_bl, atoms, label_atom, vectors)

# find the centroid of the molecule
c_x, c_y, c_z = ha.find_centroid(mol)

# translate the molecule and atoms to the centroid

for atom in atoms:
    atom.translate(-c_x, -c_y, -c_z)

# make a very big cell
mega = ha.make_mega_cell(atoms, traAN, traBN, traCN, vectors)

# get a cluster of atoms
clust = ha.make_cluster(mega, clust_rad, max_bl)

# make a list of shell atoms
shell = []
for atom in clust:
    if atom not in mol:
        shell.append(atom)


# write useful xyz
ef.write_xyz("mol.init.xyz", mol)
ef.write_xyz("clust.xyz", clust)
ef.write_xyz("shell.xyz", shell)

# EWALD
os.chdir(ewald_path)
# write inputs
ef.write_uc(name + ".uc", vectors, aN, bN, cN, atoms)
ef.write_qc(name + ".qc", clust)
ef.write_ew_in(name, "ewald.in." + name, nChk, nAt)
ef.write_seed()
# run Ewald
subprocess.call("./Ewald < ewald.in." + name, shell=True)
# read points output by Ewald
points = rf.read_points(name + ".pts-tb")
os.chdir(here)

# Make inputs
ef.write_g_temp("rl", "rl.temp", shell, [], "rl.template")
ef.write_g_temp("ml", "ml.temp", [], shell, "ml.template")
ef.write_g_temp("mh", "mh.temp", [], shell + points, "mh.template")
ef.write_g_temp("mg", "mg.temp", [], shell + points,
                "mg.template")  # only useful for CI


end_time = datetime.now()
print("ELAPSED TIME: " + str(end_time - start_time))
print("ENDING TIME: " + str(end_time))
