#!/usr/bin/env python
"""
This function uses fromage routines to find CIs without a cluster model
The function minimises the energy or the gap penalty function

"""
############################
# Developers:              #
# Dr Miguel Rivera         #
# Dr Federico Hernandez    #
############################
import numpy as np
import subprocess
import os
from datetime import datetime
from scipy.optimize import minimize

from fromage.io import read_file as rf
from fromage.utils import array_operations as ao
from fromage.utils import calc
from fromage.io.parse_config_file import bool_cast


def sequence(in_pos):
    """
    Run Gaussian calculations in parallel and write and return results
    This function is designed to work with the scipy.optimise.minimize function.
    This is why it can only receive one array of floats as input and return two
    arrays of floats. As a result some variables in this function are defined
    elsewhere in the module which is a necessary evil.

    Parameters
    ----------
    in_pos : list of floats
        Input coordinates in array form
    Returns
    -------
    en_out : float
        Combined energy or penalty function value in Hartree
    gr_out : list of floats
        Gradients of en_out in Hartree/Angstrom

    References
    ----------
    Levine, B. G., Coe, J. D. & Martinez, T. J. Optimizing conical intersections
    without derivative coupling vectors: Application to multistate
    multireference second-order perturbation theory (MS-CASPT2).
    J. Phys. Chem. B 112, 405-413 (2008).

    """
    # update geom
    with open("geom_mol.xyz", "a") as geom_m_file:
        geom_m_file.write(str(len(mol_atoms))+"\n")
        geom_m_file.write("opt"+"\n")
        for atom in ao.array2atom(mol_atoms, in_pos):
            atom_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
            geom_m_file.write(atom_str)

    # initialise calculation objects
    mh = calc.setup_calc("mh", high_level)
    mg = calc.setup_calc("mg", high_level)

    # Run the calculations as subprocesses with a maximum of 2 simultameous ones
    # at the same time. This order is optimised for the mh calculation being
    # the longest
    mh_proc = mh.run(ao.array2atom(mol_atoms, in_pos),nprocs)
    if high_level != "gaussian_cas":
        mg_proc = mg.run(ao.array2atom(mol_atoms, in_pos),nprocs)
        mg_proc.wait()
    mh_proc.wait()

    if high_level == "gaussian_cas":
        mh_en_gr = mh.read_out(in_pos)[0:3]
        mg_en_gr = (mh.read_out(in_pos)[2], mh.read_out(
            in_pos)[3], mh.read_out(in_pos)[2])
    else:
        mh_en_gr = mh.read_out(in_pos)
        mg_en_gr = mg.read_out(in_pos)

    # combine results
    en_combo = mh_en_gr[0]
    gr_combo = mh_en_gr[1]
    scf_combo = mh_en_gr[2]

    # corresponding ground state energy and gradients
    en_combo_g = mg_en_gr[0]
    gr_combo_g = mg_en_gr[1]

    # Penalty function parameters and calculation
    alpha = 0.02
    e_mean = (en_combo + en_combo_g) / 2
    e_diff = en_combo - en_combo_g
    if e_diff < 0:
        e_diff = -e_diff # account for inversion
    g_ij = e_diff**2 / (e_diff + alpha)
    en_out = e_mean + sigma * g_ij
    gr_out = 0.5 * (gr_combo + gr_combo_g) + sigma * ((e_diff**2 + 2 *
                                                       alpha * e_diff) / (e_diff + alpha)**2) * (gr_combo - gr_combo_g)

    # print some updates in the output
    out_file.write("------------------------------\n")
    global iteration
    iteration += 1
    out_file.write("Iteration: " + str(iteration) + "\n")
    out_file.write(
        "ONIOM Total energy: {:>27.8f} eV\n".format(en_combo * evconv))
    out_file.write(
        "ONIOM SCF energy: {:>29.8f} eV\n".format(scf_combo * evconv))
    out_file.write(
        "Energy grad. norm: {:>28.8f} eV/A\n".format(np.linalg.norm(gr_combo * evconv)))
    if bool_ci:
        out_file.write(
            "Penalty function value: {:>23.8f} eV\n".format(en_out * evconv))
        out_file.write("Penalty function grad. norm: {:>18.8f} eV\n".format(
            np.linalg.norm(gr_out * evconv)))

#    out_file.write("Gap: {:>42.8f} eV\n".format(
#        (en_combo - scf_combo) * evconv))
    out_file.write("Gap: {:>42.8f} eV\n".format(
        (en_combo - en_combo_g) * evconv))
    out_file.flush()
    return (en_out, gr_out)

if __name__ == '__main__':
    evconv = 27.2114  # Something in Hartree * evconv = Something in eV

    # default settings

    def_inputs = {
        "mol_file": "mol.init.xyz",
        "shell_file": "shell.xyz",
        "out_file": "fromage.out",
        "bool_ci": "1",
        "high_level": "gaussian",
        "low_level": "gaussian",
        "sigma": "3.5",
        "gtol": 1e-6,
        "nprocs": "1"}

    inputs = def_inputs.copy()

    # read user inputs
    if os.path.isfile("fromage.in"):
        new_inputs = rf.read_config("fromage.in")
        inputs.update(new_inputs)
    
    mol_file = inputs["mol_file"]
    out_file = inputs["out_file"]
    bool_ci = bool_cast(inputs["bool_ci"])
    high_level = inputs["high_level"]
    low_level = inputs["low_level"]
    gtol = inputs["gtol"]
    nprocs = inputs["nprocs"]
    # sigma is called lambda in some papers but that is a bad variable name
    # in Python
    sigma = float(inputs["sigma"])
    # output
    out_file = open(out_file, "w", 1)
    # print start time
    start_time = datetime.now()
    out_file.write("STARTING TIME: " + str(start_time) + "\n")

    iteration = 0

    # clean up the last output
    if os.path.exists("geom_mol.xyz"):
        subprocess.call("rm geom_mol.xyz", shell=True)
    if os.path.exists("geom_cluster.xyz"):
        subprocess.call("rm geom_cluster.xyz", shell=True)

    # read initial coordniates
    mol_atoms = rf.read_xyz(mol_file)[0]

    # make the initial coordinates into a flat list
    atoms_array = []
    for atom in mol_atoms:
        atoms_array.append(atom.x)
        atoms_array.append(atom.y)
        atoms_array.append(atom.z)

    # make the list into an array
    atoms_array = np.array(atoms_array)
    res = minimize(sequence, atoms_array, jac=True,
                   options={'disp': True, 'gtol': gtol})

    out_file.write("DONE\n")
    end_time = datetime.now()
    out_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    out_file.write("ENDING TIME: " + str(end_time) + "\n")
    out_file.close()
