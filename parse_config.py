#!/usr/bin/env python
import read_file as rf

name = "config"

cell_file = name + "xyz"
# name of the cp2k file with population information
cp2k_file = False
# name of the Gaussian log file with population information
mol_pop_file = False
# kind of population in the Gaussian file 0:Mulliken 1:ESP
mol_pop_kind = 0
# maximum bond length when defining a molecule
max_bl = 1.7
# label of an atom which will be part of the quantum cluster
# warning: [0,N-1], not [1,N]
label_atom = 0
# the number of checkpoints in region 1
nChk = 1000
# the number of constrained charge atoms
# i.e. atoms in regions 1 and 2
nAt = 500
# Ewald will multiply the unit cell in the direction
# of the a, b or c vector 2N times (N positive and N negative)
aN = 2
bN = 2
cN = 2
# Population analysis method if pertinent
# Mulliken(0) Hirshfeld(1) RESP(2)
cp2k_pop_method = 0
# the cluster will be of all molecules with atoms less than
# clust_rad away from the centre of the central molecule
clust_rad = 5
# how many times the input cluster needs to be repeated along each vector
# positively and negatively to be able to contain the cluster to select.
# the supercluster ends up being (1+2*traAN)*(1+2*traBN)*(1+2*traCN) times
# bigger
traAN = 2
traBN = 2
traCN = 2
# Self Consistent Ewald Gaussian template
sc_temp = False
# Self Consistent Ewald atom kind 0:Mulliken 1:ESP
sc_kind = 0
# Self Consistent Ewald deviation tolerance
dev_tol = 0.001
# Ewald embedding, use 0 for false
ewe = True
# Cluster self consistent template for mol
csc_temp_h = False
# Cluster self consistent template for shell
csc_temp_l = ""

dic = rf.read_config(name)

print(dic)
