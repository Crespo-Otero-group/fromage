#!/usr/bin/env python
"""
Minimises the energy or the gap penalty function of a molecule in a cluster.
Also performs Surface Hopping dynamics at CASSCF level with OpenMolcas

Template files are needed as well as an xyz file for the molecule and another
one for the surrounding molecules. Overall the use of subprocess is ugly as it
is repeated 3 or 4 times but it was found to handle memory better than Pool
when interfacing with Gaussian.

Energy units are Hartree inside the program but are printed in eV. Distances are
kept as Angstrom throughout and are converted from Bohr if necessary before
reaching this module

"""
import numpy as np
import subprocess
import os
from datetime import datetime
from scipy.optimize import minimize

from fromage.io import read_file as rf
from fromage.utils import array_operations as ao
from fromage.utils import calc
from fromage.utils import fro_dyn as fd
from fromage.utils import vib_analysis as va
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
    # initialise calculation objects
    rl = calc.setup_calc("rl", low_level)
    ml = calc.setup_calc("ml", low_level)
    mh = calc.setup_calc("mh", high_level)
    if bool_ci:
        mg = calc.setup_calc("mg", high_level)

    # Run the calculations as subprocesses with a maximum of 2 simultameous ones
    # at the same time. This order is optimised for the mh calculation being
    # the longest

    all_pos = np.concatenate((in_pos, fixed_atoms_array), axis = 0)

    if low_level == "fomo-ci" or low_level == "mopac" and at_reparam is not None:
        rl_proc = rl.run(ao.array2atom(all_atoms, all_pos),nprocs,at_reparam)
        rl_proc.wait()
    else:
        rl_proc = rl.run(atoms = ao.array2atom(all_atoms, all_pos), nprocs = nprocs)
        rl_proc.wait()

    # Get the charges and use them for the model region

    rl_charges_array = rl.read_charges()

    if high_level == "fomo-ci" or high_level == "mopac" and at_reparam is not None:
        mh_proc = mh.run(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                         ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),
                         nprocs, at_reparam)
        if bool_ci:
            mg_proc = mg.run(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                      ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),
                      nprocs,at_reparam)
    else:
        mh_proc = mh.run(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                         ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),nprocs)
        mh_proc.wait()
        if bool_ci and high_level != "gaussian_cas":
            mg_proc = mg.run(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                             ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),nprocs)
            mg_proc.wait()
#    mh_proc.wait()
    if low_level == "fomo-ci" or low_level == "mopac" and at_reparam is not None:
        ml_proc = ml.run(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                         ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),
                         nprocs,at_reparam)
        ml_proc.wait()
    else:
        ml_proc = ml.run(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                         ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]), nprocs)
        ml_proc.wait()

    # read results. Each x_en_gr is a tuple (energy,gradients,scf_energy)
    rl_en_gr = rl.read_out(in_pos, 
                           in_mol = all_flex_atoms, # FJH here I giving QM and flex QM' region instead of mol_atoms
                           in_shell = fixed_atoms, # FJH here I giving only the fixed atoms instead of all the shell_atoms
                           natoms_flex = natoms_flex)
    ml_en_gr = ml.read_out(in_pos[:dim_qm], natoms_flex = natoms_flex)
#
    if high_level == "gaussian_cas":
        mh_en_gr = mh.read_out(in_pos[:dim_qm],natoms_flex = natoms_flex)[0:3]
        if bool_ci:
            mg_en_gr = (mh.read_out(in_pos[:dim_qm],natoms_flex = natoms_flex)[2], mh.read_out(
                in_pos[:dim_qm], natoms_flex = natoms_flex)[3], mh.read_out(
                    in_pos[:dim_qm], natoms_flex = natoms_flex)[2])
    else:
        mh_en_gr = mh.read_out(in_pos[:dim_qm], natoms_flex = natoms_flex)
        if bool_ci:
            mg_en_gr = mg.read_out(in_pos[:dim_qm], natoms_flex = natoms_flex)

    # combine results
    en_combo = rl_en_gr[0] - ml_en_gr[0] + mh_en_gr[0]
    scf_combo = rl_en_gr[2] - ml_en_gr[2] + mh_en_gr[2]

    if single_point:
        gr_combo = 0
    else:
        gr_combo = rl_en_gr[1] - ml_en_gr[1] + mh_en_gr[1]

    # if linker atoms are included, gr_combo has to be defined as:
      # where J is the Jacobian that can be easily defined according to
      # the Morokuma's definition https://doi.org/10.1021/cr5004419
      #gr_combo = rl_en_gr[1] - ml_en_gr[1] x Jac + mh_en_gr[1] x Jac

    if bool_ci:
        # corresponding ground state energy and gradients
        en_combo_g = rl_en_gr[0] - ml_en_gr[0] + mg_en_gr[0]
        gr_combo_g = rl_en_gr[1] - ml_en_gr[1] + mg_en_gr[1]
        #gr_combo = rl_en_gr[1] - ml_en_gr[1] x Jac + mh_en_gr[1] x Jac

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
        e_diff = 0

    global iteration
    iteration += 1
    _write_calc_info(out_file = out_file,
                     mh_en_gr = mh_en_gr,
                     ml_en_gr = ml_en_gr,
                     rl_en_gr = rl_en_gr,
                     en_combo = en_combo,
                     gr_combo = gr_combo,
                     scf_combo = scf_combo,
                     evconv = evconv,
                     iteration = iteration,
                     en_out = en_out,
                     gr_out = gr_out,
                     e_diff = e_diff,
                     bool_ci = bool_ci)

    return (en_out, gr_out)

#start_trajectory(dyn_array,inputs,mol_atoms,flex_atoms,fixed_atoms)

def start_trajectory(geometry,dyn_sett,mol_atoms,shell_atoms):
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
    traj = fd.Trajectory(in_params, mol_atoms, shell_atoms)

    traj.run_dynamics()

    return None

###########################################################################
################################### FJH ###################################
def set_newtonx(atoms_array,inputs):
    """
    This subroutine prepare all the environments and files to use fromage
    as a third-party program of Newton-X for the calculation of spectra
    and dynamics
    """
    natoms, nstates, state = nx.read_nx_control()
    nx.newtonx_sequence(atoms_array,inputs,natoms,nstates,state)

    return None
############################################################################

def start_normal_modes(inputs,mol_atoms,QM_natoms,natoms_flex,flex_atoms,fixed_atoms,fixed_atoms_array):
    """
    """
    out_file.write("A calculation of ONIOM normal modes has been requested\n")
    if os.path.exists("geom_mol.xyz"):
        geom_mol = True
        out_file.write("geom_mol.xyz file has been found\n")
        out_file.write("The cluster geometry will be set up from the last geometry in geom_mol.xyz \n")
    else:
        out_file.write("The cluster geometry will be set up from the geometries in mol.init.xyz and shell_flex.xyz \n")
    at_symbols, at_pos = va.get_xyz_cluster(QM_natoms, natoms_flex, flex_atoms, geom_mol)

    in_params = va.initNmodesParams(at_symbols,at_pos,inputs)

    Nmodes = va.NormalModes(in_params,natoms_flex,mol_atoms,flex_atoms,fixed_atoms,fixed_atoms_array)

    # Compute the ONIOM normal modes
    Nmodes.compute_nmodes()

    return None

def _write_head(out_file):
    """
    """
    # print start time
    start_time = datetime.now()
    out_file.write("STARTING TIME: " + str(start_time) + "\n")
    out_file.write("" "\n")
    out_file.write("************************************************" "\n")
    out_file.write(" Find the bug between the code and the output " "\n")
    out_file.write("\n")
    out_file.write("If you see something that it doesn't look right" "\n")
    out_file.write("          See it, Say it, Sort it...           " "\n")
    out_file.write("\n")
    out_file.write("************************************************" "\n")
    return start_time

def _write_calc_info(out_file,
                     mh_en_gr,
                     ml_en_gr,
                     rl_en_gr,
                     en_combo,
                     gr_combo,
                     scf_combo,
                     evconv,
                     iteration,
                     en_out,
                     gr_out,
                     e_diff,
                     bool_ci = None):
    """
    print some updates in the output
    """
    out_file.write("------------------------------\n")
    out_file.write("Iteration: " + str(iteration) + "\n")
    out_file.write("Real low energy: {:>30.8f} eV\n".format(
        rl_en_gr[0] * evconv))
    out_file.write("Model low energy: {:>29.8f} eV\n".format(
        ml_en_gr[0] * evconv))
    out_file.write("Model high energy: {:>28.8f} eV\n".format(
        mh_en_gr[0] * evconv))
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
        out_file.write("Gap: {:>42.8f} eV\n".format(
            e_diff*evconv))
    else:
        out_file.write("Gap: {:>42.8f} eV\n".format(
            (en_combo - scf_combo) * evconv))
        out_file.flush()

    return

def _write_tail(start_time,out_file):
    """
    Writes the time info when the optimization process
    or dynamics is finished
    """
    out_file.write("DONE\n")
    end_time = datetime.now()
    out_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    out_file.write("ENDING TIME: " + str(end_time) + "\n")
    out_file.close()

    return None


if __name__ == '__main__':

    evconv = 27.2114  # Something in Hartree * evconv = Something in eV

    # default settings

    def_inputs = {
        "mol_file": "mol.init.xyz",
        "shell_file_fixed": "shell_fixed.xyz",
        "shell_file_flex": "shell_flex.xyz",
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
        "relax": "0",
        "at_reparam": "0",
        "natoms_flex": "0",
        "normal_modes" : "0"} 

    inputs = def_inputs.copy()

    # read user inputs
    if os.path.isfile("fromage.in"):
        new_inputs = rf.read_config("fromage.in")
        inputs.update(new_inputs)

    out_file = inputs["out_file"]
#
    # output
    out_file = open(out_file, "w", 1)
    # write head in the output file
    start_time = _write_head(out_file)
#
    natoms_flex = int(inputs["natoms_flex"])
    relax = bool_cast(inputs["relax"])
    normal_modes = bool_cast(inputs["normal_modes"])
    if natoms_flex > 0 and not relax and not normal_modes:
        out_file.write(" "+ "\n")
        out_file.write("You are trying to run fromage with the relaxation of the QMprime region option OFF"+ "\n")
        out_file.write("Include *relax* keyword in fromage.in"+ "\n")
        out_file.write("fromage is dying now :-( "+ "\n")
        import sys
        sys.exit()
    elif (relax and natoms_flex == 0) or (normal_modes and natoms_flex == 0) :
        out_file.write(" "+ "\n")
        out_file.write("You are trying to run either fromage with the relaxation of the QMprime region option ON"+ "\n")
        out_file.write("You are trying to run either fromage with ONIOM normal modes"+ "\n")
        out_file.write("However, the number of flexible atoms set in fromage.in is %s" % (natoms_flex) + "\n")
        out_file.write(" Include the keyword *natoms_flex* along with the number of flexible atoms in fromage.in and re run"+ "\n")
        out_file.write("fromage is dying now :-( "+ "\n")
        import sys
        sys.exit()

    mol_file = inputs["mol_file"]
    shell_file_flex = inputs["shell_file_flex"]
    shell_file_fixed = inputs["shell_file_fixed"]
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
    # If so, the atom number is collected and a "w" symbol is added next to 
    # the atom symbol in the coordinates added to the FOMO-CI input.
    if "at_reparam" in inputs.keys():
         at_reparam = [] 
         at_reparam = [int(x) for x in inputs["at_reparam"]]
         at_reparam = np.array(at_reparam)
    else:
        at_reparam = None

    if nprocs=="1":
        out_file.write("If Q-Chem, Molcas, NWChem or MOPAC are to be used, have in mind that" "\n")
        out_file.write("the default number of cores are asked for the calculation,")
        out_file.write("regardless what you have asked in your submission script file: " + "nprocs=" + str(nprocs) + "\n")
        out_file.write("" "\n")

    iteration = 0

    # clean up the last output
    if normal_modes is None:
        if os.path.exists("geom_mol.xyz"):
            subprocess.call("rm geom_mol.xyz", shell=True)
        if os.path.exists("geom_cluster.xyz"):
            subprocess.call("rm geom_cluster.xyz", shell=True)

    # read initial coordniates
    mol_atoms = rf.read_xyz(mol_file)[0]   
    # read atoms for the flexible shell
    flex_atoms = rf.read_xyz(shell_file_flex)[0]
    # read atoms for the fixed shell
    fixed_atoms = rf.read_xyz(shell_file_fixed)[0]
    # Define the entire flexible region
    all_flex_atoms = mol_atoms + flex_atoms
    all_atoms = mol_atoms + flex_atoms + fixed_atoms
    # make the initial coordinates into a flat list
    atoms_array = []
    fixed_atoms_array = []
    dyn_array = []
    for atom in mol_atoms:
        dyn_array.append([ Element(atom.at_num).getSymbol(), atom.x, atom.y, atom.z ])
        atoms_array.append(atom.x)
        atoms_array.append(atom.y)
        atoms_array.append(atom.z)
    for atom in flex_atoms:
        dyn_array.append([ Element(atom.at_num).getSymbol(), atom.x, atom.y, atom.z ])
        atoms_array.append(atom.x)
        atoms_array.append(atom.y)
        atoms_array.append(atom.z)
    for atom in fixed_atoms:
        fixed_atoms_array.append(atom.x)
        fixed_atoms_array.append(atom.y)
        fixed_atoms_array.append(atom.z)
#        charges_array.append(atom.q)	#Initial charges
    # make the list into an array
    atoms_array = np.array(atoms_array)
    fixed_atoms_array = np.array(fixed_atoms_array)
#    charges_array = np.array(charges_array)
    dim_flex = int(3*natoms_flex)
    QM_natoms = len(mol_atoms)
    dim_qm = int(3*QM_natoms)
    
    if single_point:
        out_file.write("A single point calculation has been requested\n")
        sequence(atoms_array)
    elif dynamics:
        out_file.write("A dynamics calculation has been requested\n")
        res = start_trajectory(dyn_array,inputs,mol_atoms,shell_atoms)
#        res = start_trajectory(dyn_array,inputs,mol_atoms,flex_atoms,fixed_atoms) # FJH
    elif relax:
        out_file.write("An ONIOM optimization has been requested\n")
        res = minimize(sequence, atoms_array, jac=True,
                       options={'disp': True, 'gtol': gtol})
    if normal_modes:
        res = start_normal_modes(inputs,mol_atoms,QM_natoms,natoms_flex,flex_atoms,fixed_atoms,fixed_atoms_array)

    _write_tail(start_time,out_file)

