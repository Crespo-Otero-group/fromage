"""Minimises the energy or the gap penalty function of a molecule in a cluster.

Template files are needed as well as an xyz file for the molecule and another
one for the surrounding molecules. Overall the use of subprocess is ugly as it
is repeated 3 or 4 times but it was found to handle memory better than Pool
when interfacing with Gaussian.

"""
import numpy as np
import subprocess
import os
import read_file as rf
import edit_file as ef
from atom import Atom
from scipy.optimize import minimize
from datetime import datetime


# output
out_file = open("cryspy.out", "w", 1)
# print start time
start_time = datetime.now()
out_file.write("STARTING TIME: " + str(start_time)+"\n")


def array2atom(template, pos):
    """
    Turn an array of the form x1, y1, z1, x2, y2, z2 etc. into a list of Atom
    objects

    Parameters
    ----------
    template : list of Atom objects
        A list of the same length of the desired one used to determine the
        elements of the atoms
    pos : list of floats
        List of coordinates of the form x1, y1, z1, x2, y2, z2 etc.
    Returns
    ----------
    out_atoms : list of Atom objects
        Resulting atoms

    """
    sliced_pos = [pos[i:i + 3] for i in range(0, len(pos), 3)]
    out_atoms = []
    for atom in zip(template, sliced_pos):
        new_atom = Atom(atom[0].elem, atom[1][0], atom[1][1], atom[1][2], 0)
        out_atoms.append(new_atom)
    return out_atoms


def treat_chk(positions, name, print_bool, in_mol, in_shell):
    """
    Analyse a Gaussian .chk file while printing geometry updates

    Parameters
    ----------
    positions : list of floats
        List of atomic coordinates
    name : string
        Name of the .chk file to treat without the extension the extension
    print_bool : bool
        If true, update the geometry files with the one found in .chk
    in_mol : list of Atom objects
        Atoms in the inner region
    in_shell : list of Atom objects
        Atoms in the middle region
    Returns
    ----------
    energy : float
        Energy calculated by Gaussian in Hartree
    gradients : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
    scf_energy : float
        The ground state energy in Hargree

    """
    # stdout=FNULL to not have to read the output of formchk
    FNULL = open(os.devnull, 'w')
    subprocess.call("formchk " + name + ".chk", stdout=FNULL, shell=True)
    energy, gradients_b, scf_energy = rf.read_fchk(name + ".fchk")
    # fix gradients units
    gradients = gradients_b * bohrconv
    # update the geometry log
    if print_bool == True:
        # only the inner region
        with open("geom_mol.xyz", "a") as geom_m_file:
            geom_m_file.write(str(len(in_mol)) + "\n")
            geom_m_file.write("E_" + name + "= " + str(energy) + "\n")
            for atom in array2atom(in_mol, positions):
                atomStr = "{:>6} {:10.9f} {:10.6f} {:10.9f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_m_file.write(atomStr)
        # the inner and middle regions
        with open("geom_cluster.xyz", "a") as geom_c_file:
            geom_c_file.write(str((len(positions) / 3) + len(in_shell)) + "\n")
            geom_c_file.write("E_" + name + "= " + str(energy) + "\n")
            for atom in array2atom(in_mol, positions):
                atomStr = "{:>6} {:10.9f} {:10.6f} {:10.9f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_c_file.write(atomStr)
            for atom in in_shell:
                atomStr = "{:>6} {:10.9f} {:10.9f} {:10.9f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_c_file.write(atomStr)

    # truncate gradients if too long
    gradients = gradients[:len(in_mol * 3)]

    return (energy, gradients, scf_energy)


def sequence(in_pos):
    """
    Run Gaussian calculations in parallel and write and return results

    This function is designed to work with the scipy.optimise.minimize function.
    This is why it can only receive one array of floats as input and return two
    arrays of floats. As a result some variables in this function are defined
    elsewhere in the module which is a necessary evil.

    Parameter
    ----------
    in_pos : list of floats
        Input coordinates in array form
    Returns
    ----------
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
    # write Gaussian .com files for all calculations
    ef.write_gauss("rl", "rl.com", array2atom(
        mol_atoms, in_pos), [], "rl.temp")
    ef.write_gauss("ml", "ml.com", array2atom(
        mol_atoms, in_pos), [], "ml.temp")
    ef.write_gauss("mh", "mh.com", array2atom(
        mol_atoms, in_pos), [], "mh.temp")
    if bool_ci:
        ef.write_gauss("mg", "mg.com", array2atom(
            mol_atoms, in_pos), [], "mg.temp")

    # Run the calculations as subprocesses with a maximum of 2 simultameous ones
    # at the same time. This order is optimised for the mh calculation being
    # the longest
    p_mh = subprocess.Popen("g09 mh.com", shell=True)
    p_rl = subprocess.Popen("g09 rl.com", shell=True)
    p_rl.wait()
    p_ml = subprocess.Popen("g09 ml.com", shell=True)
    p_ml.wait()
    if bool_ci:
        p_mg = subprocess.Popen("g09 mg.com", shell=True)
        p_mg.wait()
    p_mh.wait()

    # read results. Each x_en_gr is a tuple (energy,gradients,scf_energy)
    rl_en_gr = treat_chk(in_pos, "rl", True, mol_atoms, shell_atoms)
    ml_en_gr = treat_chk(in_pos, "ml", False, mol_atoms, shell_atoms)
    mh_en_gr = treat_chk(in_pos, "mh", False, mol_atoms, shell_atoms)
    if bool_ci:
        mg_en_gr = treat_chk(in_pos, "mg", False, mol_atoms, shell_atoms)

    # combine results
    en_combo = rl_en_gr[0] - ml_en_gr[0] + mh_en_gr[0]
    gr_combo = rl_en_gr[1] - ml_en_gr[1] + mh_en_gr[1]
    scf_combo = rl_en_gr[2] - ml_en_gr[2] + mh_en_gr[2]

    if bool_ci:
        # corresponding ground state energy and gradients
        en_combo_g = rl_en_gr[0] - ml_en_gr[0] + mg_en_gr[0]
        gr_combo_g = rl_en_gr[1] - ml_en_gr[1] + mg_en_gr[1]

        # Penalty function parameters and calculation
        alpha = 0.02
        # sigma is called lambda in some papers but that is a bad variable name
        # in Python
        sigma = 3.5
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
    out_file.write("Iteration: " + str(iteration)+"\n")
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
    if bool_ci == True:
        out_file.write("Penalty function value: {:>23.8f} eV\n".format(en_combo * evconv))
        out_file.write("Penalty function grad. norm: {:>18.8f} eV\n".format(
            np.linalg.norm(gr_combo * evconv)))
    out_file.write("Gap: {:>43.8f} eV\n".format(
        (en_combo - scf_combo) * evconv))
    out_file.flush()
    return (en_out, gr_out)

evconv = 27.2114  # Something in Hartree * evcomv = Something in eV
bohrconv = 1.88973  # Something in Angstrom 8 bohrconv = Something in Bohr

iteration = 0

# is this a CI calculation?
bool_ci = False

# clean up the last output
if os.path.exists("geom_mol.xyz"):
    subprocess.call("rm geom_mol.xyz", shell=True)
if os.path.exists("geom_cluster.xyz"):
    subprocess.call("rm geom_cluster.xyz", shell=True)


# read initial coordniates
mol_atoms = rf.read_xyz("mol.init.xyz")[0]

# read shell atoms
shell_atoms = rf.read_xyz("shell.xyz")[0]

# make the initial coordinates into a flat list
atoms_array = []
for atom in mol_atoms:
    atoms_array.append(atom.x)
    atoms_array.append(atom.y)
    atoms_array.append(atom.z)

# make the list into an array
atoms_array = np.array(atoms_array)


res = minimize(sequence, atoms_array, jac=True,
               options={'disp': True})

out_file.write("DONE\n")
end_time = datetime.now()
out_file.write("ELAPSED TIME: " + str(end_time - start_time)+"\n")
out_file.write("ENDING TIME: " + str(end_time)+"\n")
out_file.close()
