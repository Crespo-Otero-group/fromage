#!/usr/bin/env python
"""
Fromage function to perform MECI search with the penalty function method
or Surface Hopping dynamics in the gas phase

Template files are needed as well as an xyz file for the molecule.

Energy units are Hartree inside the program but are printed in eV. Distances are
kept as Angstrom throughout and are converted from Bohr if necessary before
reaching this module

"""
##########################
#       Written by       #
# Dr Federico Hernandez  #
##########################
import numpy as np
import subprocess
import os
from datetime import datetime
from scipy.optimize import minimize

from fromage.io import read_file as rf
from fromage.utils import array_operations as ao
from fromage.utils import calc
from fromage.gas_phase import fro_gp_dyn as fd
from fromage.io.parse_config_file import bool_cast
from fromage.dynamics.periodic_table import Element


def sequence(in_pos):
    """
    Run Electronic Structure calculations in parallel and write and return results
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

    Updated by Federico Hernandez 27-10-2022
    """
    # update geom
    with open("geom_mol.xyz", "a") as geom_m_file:
        geom_m_file.write(str(len(mol_atoms))+"\n")
        geom_m_file.write("opt MECI"+"\n")
        for atom in ao.array2atom(mol_atoms, in_pos):
            atom_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
            geom_m_file.write(atom_str)

    # initialise calculation objects
    mh = calc.setup_calc("mh", high_level)
    if bool_ci:
        mg = calc.setup_calc("mg", high_level)

    # Run the calculations as subprocesses with a maximum of 2 simultameous ones
    # at the same time. This order is optimised for the mh calculation being
    # the longest

    if high_level == "fomo-ci" or high_level == "mopac" and at_reparam is not None:
        mh_proc = mh.run(atoms = ao.array2atom(mol_atoms, in_pos), 
                         nprocs = nprocs, at_reparam = at_reparam)
        if bool_ci:
            mg_proc = mg.run(atoms = ao.array2atom(mol_atoms, in_pos), 
                             nprocs = nprocs, at_reparam = at_reparam)
    else:
        mh_proc = mh.run(atoms = ao.array2atom(mol_atoms, in_pos), nprocs = nprocs)
        mh_proc.wait()
        if bool_ci and high_level != "gaussian_cas":
            mg_proc = mg.run(atoms = ao.array2atom(mol_atoms, in_pos), nprocs = nprocs)
            mg_proc.wait()
#        mh_proc.wait()
#
    if high_level == "gaussian_cas":
        mh_en_gr = mh.read_out(in_pos)[0:3]
        if bool_ci:
            mg_en_gr = (mh.read_out(in_pos)[2], mh.read_out(
                in_pos)[3], mh.read_out(in_pos)[2])
    else:
        mh_en_gr = mh.read_out(in_pos)
        if bool_ci:
            mg_en_gr = mg.read_out(in_pos)

    # combine results
    en_combo = mh_en_gr[0]
    gr_combo = mh_en_gr[1]
    scf_combo = mh_en_gr[2]

    if bool_ci:
        # corresponding ground state energy and gradients
        en_combo_g = mg_en_gr[0]
        gr_combo_g = mg_en_gr[1]

        # Penalty function parameters and calculation
        alpha = 0.02
        e_mean = (en_combo + en_combo_g) / 2
        e_diff = en_combo - en_combo_g
        g_ij = e_diff**2 / (e_diff + alpha)
        en_out = e_mean + sigma * g_ij
        gr_out = 0.5 * (gr_combo + gr_combo_g) + sigma * ((e_diff**2 + 2 *
                                                           alpha * e_diff) / (e_diff + alpha)**2) * (gr_combo - gr_combo_g)
    else:
        en_out = en_combo
        gr_out = gr_combo

    # print some updates in the output
    out_file.write("------------------------------\n")
    global iteration
    iteration += 1
    out_file.write("Iteration: " + str(iteration) + "\n")
    out_file.write(
        "Excited state energy: {:>27.8f} eV\n".format(en_combo * evconv))
    out_file.write(
        "Ground state energy: {:>29.8f} eV\n".format(scf_combo * evconv))
    out_file.write(
        "Energy grad. norm: {:>28.8f} eV/A\n".format(np.linalg.norm(gr_combo * evconv)))
    if bool_ci:
        out_file.write(
            "Penalty function value: {:>23.8f} eV\n".format(en_out * evconv))
        out_file.write("Penalty function grad. norm: {:>18.8f} eV\n".format(
            np.linalg.norm(gr_out * evconv)))
        out_file.write("Gap: {:>42.8f} eV\n".format(
            (en_combo - en_combo_g) * evconv))
    else:
        out_file.write("Gap: {:>42.8f} eV\n".format(
            (en_combo - scf_combo) * evconv))
        out_file.flush()
    return (en_out, gr_out)

def start_trajectory(geometry, dyn_sett, mol_atoms):
    # Read initial velocities from file
    in_vel = rf.read_velocities(dyn_sett["vel_file"])
    # Read the gradient of the step previous the dynamics crashed
    if dyn_restart:
        curr_step, Eini, prev_grad = rf.read_dyn_restart("dyn_restart")

    # Create Trajectory object with dynamics info and initial conditions
    atomic_symbols = [ x[0] for x in geometry ]
    in_pos = [ [x[1], x[2], x[3]] for x in geometry ]
    if dyn_restart:
        in_params = fd.initTrajParams(atomic_symbols, in_pos, in_vel, dyn_sett, curr_step, Eini, prev_grad)
    else:
        in_params = fd.initTrajParams(atomic_symbols, in_pos, in_vel, dyn_sett)
    traj = fd.Trajectory(in_params, mol_atoms)

    traj.run_dynamics()

    return None


if __name__ == '__main__':

    evconv = 27.2114  # Something in Hartree * evconv = Something in eV

    # default settings

    def_inputs = {
        "mol_file": "mol.init.xyz",
        "shell_file": "shell.xyz",
        "out_file": "fromage.out",
        "bool_ci": "0",
        "high_level": "gaussian",
        "low_level": "gaussian",
        "nprocs": "1",
        "sigma": "3.5",
        "gtol": "1e-5",
        "single_point": "0",
        "dynamics": "0",
        "dyn_restart": "0",
        "relax_qmprime": "0",
        "at_reparam": "0",
        "natoms_flex": "0"} 

    inputs = def_inputs.copy()

    # read user inputs
    if os.path.isfile("fromage.in"):
        new_inputs = rf.read_config("fromage.in")
        inputs.update(new_inputs)

########### FJH #################################
    out_file = inputs["out_file"]
#
    # output
    out_file = open(out_file, "w", 1)
    # print start time
    start_time = datetime.now()
    out_file.write("STARTING TIME: " + str(start_time) + "\n")
#
    natoms_flex = bool_cast(inputs["natoms_flex"])
    relax_qmprime = bool_cast(inputs["relax_qmprime"])
    if relax_qmprime:
        out_file.write(" "+ "\n")
        out_file.write("You are trying to run a gas phase calculation with the relaxation of the QMprime region ON"+ "\n")
        out_file.write("You should run *fro_run_flex.py* instead"+ "\n")
        out_file.write("fromage is dying now :-("+ "\n")
        import sys 
        sys.exit()
#################################################
    mol_file = inputs["mol_file"]
    bool_ci = bool_cast(inputs["bool_ci"])
    high_level = inputs["high_level"]
    low_level = inputs["low_level"]
    nprocs = inputs["nprocs"]
    gtol = float(inputs["gtol"])
    single_point = bool_cast(inputs["single_point"])
    dynamics = bool_cast(inputs["dynamics"])
    dyn_restart = bool_cast(inputs["dyn_restart"])
    # sigma is called lambda in some papers but that is a bad variable name
    # in Python
    sigma = float(inputs["sigma"])
    # Check if the are are atoms to be reparametrised for a FOMO-CI calc. 
    #If so, the atom number is collected and a "w" symbol is added next to 
    # the atom symbol in the coordinates added to the FOMO-CI input.
    if "at_reparam" in inputs.keys():
         at_reparam = [] 
         at_reparam = [int(x) for x in inputs["at_reparam"]]
         at_reparam = np.array(at_reparam)
    else:
        at_reparam = None

    if nprocs=="1":
        out_file.write("If Q-Chem, Molcas, NWChem or MOPAC are used, have in mind that" "\n")
        out_file.write("the default number of cores are asked for the calculation,")
        out_file.write("regardless what you have asked in your submission script file: " + "nprocs=" + str(nprocs) + "\n")
        out_file.write("" "\n")

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
    dyn_array = []
    for atom in mol_atoms:
        dyn_array.append([ Element(atom.at_num).getSymbol(), atom.x, atom.y, atom.z ])
        atoms_array.append(atom.x)
        atoms_array.append(atom.y)
        atoms_array.append(atom.z)
#
    # make the list into an array
    atoms_array = np.array(atoms_array)
    if single_point:
        sequence(atoms_array)
    elif dynamics:
        res = start_trajectory(dyn_array, inputs, mol_atoms)
    else:
        res = minimize(sequence, atoms_array, jac=True,
                       options={'disp': True, 'gtol': gtol})

    out_file.write("DONE\n")
    end_time = datetime.now()
    out_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    out_file.write("ENDING TIME: " + str(end_time) + "\n")
    out_file.close()
