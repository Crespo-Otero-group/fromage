#!/usr/bin/env python
"""This module is a tool to prepare template and input files for cryspy.

You will need:
- config file
- .xyz file of a unit cell
- cp2k or gaussian output files for a population analysis at high level and low
level
- .template files for ml, mh, mg, rl calculations

And receive
- .temp files for ml, mh, mg, rl
- mol.init.xyz containing the selected molecule(s)
- clust.xyz containing the full molecular cluster
- shell.xyz containing the cluster without the selected molecule(s)
"""
import subprocess
import os
import numpy as np
from datetime import datetime
from copy import copy

from cryspy.io import read_file as rf
from cryspy.io import edit_file as ef
from cryspy.io import parse_config_file as pcf
from cryspy.utils import handle_atoms as ha
from cryspy.scripts.assign_charges import assign_charges
from cryspy.utils.mol import Mol


def run_ewald(in_name, in_mol, in_atoms, in_vectors, in_nAt=500, in_aN=2, in_bN=2, in_cN=2, in_nChk=1000):
    """
    Perform an Ewald calculation in the working directory

    Parameters
    ----------
    in_name : str
        Name of the calculation
    in_mol : list of Atom objects
        Atoms which constitute the quantum cluster
    in_atoms : list of Atom objects
        Atoms in the unit cell
    in_vectors : 3x3 numpy matrix
        Lattice vectors
    in_nAt : int
        Amount of atoms with fixed charge
    in_aN, in_bN, in_cN : ints
        Amount of multiplications in postitive and negative directions of the unit cell
    in_nChk : int
        Number of random sampling points in the qc
    """
    FNULL = open(os.devnull, 'w')
    ef.write_uc(in_name + ".uc", in_vectors, in_aN, in_bN, in_cN, in_atoms)
    ef.write_qc(in_name + ".qc", in_mol)
    ef.write_ew_in(in_name, "ewald.in." + in_name, in_nChk, in_nAt)
    ef.write_seed()
    # run Ewald
    subprocess.call("Ewald < ewald.in." + in_name, stdout=FNULL, shell=True)
    return


if __name__ == '__main__':

    # a few functions to only be used in main
    def populate_cell(in_mol, program, pop_file, method):
        """
        Assign charge to the atoms from a unit cell

        Don't import this.
        If you want to assign charges go straign to assign_charges.py and use
        those utilities.

        Parameters
        ----------
        in_mol : Mol object
            Make sure these atoms form a unit cell.
        program : str
            Make it "gaussian" or "cp2k"
        pop_file : str
            The corresponding file to read charges from
        method : str
            Acceptable strings are "esp", "mulliken" and "hirshfeld"

        """
        if program.lower() == "cp2k":
            charges = rf.read_cp2k(pop_file, method)[0]
            output_file.write("Read " + str(len(in_mol)) +
                              " charges in cp2k_file\n")
            # in case there are more charges than atoms
            charges = charges[:len(in_mol)]
            # correct charges if they are not perfectly neutral
            if sum(charges) != 0.0:
                output_file.write("Charge correction: " +
                                  str(sum(charges)) + "\n")
                charges[-1] -= sum(charges)

        if program.lower() == "gaussian":
            mol_char = rf.mol_from_gauss(pop_file, pop=method)

            charges = [i.q for i in mol_char]
            # correct charges if they are not perfectly neutral
            if sum(charges) != 0.0:
                output_file.write("Charge correction: " +
                                  str(sum(charges)) + "\n")
                mol_char[-1].q -= sum(mol_char)

            # assign charges to the rest of the cell
            assign_charges(mol_char, in_mol)
        return

    output_file = open("prep.out", "w")

    def ec_loop(in_shell, in_mol, damping=0):
        """Do a single iteration Gaussian loop with no Ewald

        Beware! This function changes the charge values of in_atoms and by
        extension in_mol. The point is to run this until convergence.

        Parameters:
        -----------
        in_shell : list of Atom objects
            They should form a cell
        in_mol : list of Atom objects
            A list of atoms which are a part of in_atoms and form the quantum
            cluster
        damping : float
            Underrelaxation damping constant
        Returns:
        --------
        deviation : float
            Average deviation of charges from in_mol after an ewald calculation
            followed by a Gaussian calculation

        """
        global sc_loop
        sc_name = "sc_" + name
        sc_loop += 1

        # Initial charges in mol
        old_charges = [atom.q for atom in mol]

        # Calculate new charges
        ef.write_gauss(sc_name + ".com", mol, in_shell,
                       sc_temp, proj_name=sc_name)
        subprocess.call("g09 " + sc_name + ".com", shell=True)

        intact_charges, new_energy, char_self, char_int = rf.read_g_char(
            sc_name + ".log", high_pop_method, debug=True)
        # Correct charges if they are not perfectly neutral
        if sum(intact_charges) != 0.0:
            # intact_charges[-1] -= sum(new_charges)
            temp_correct = sum(intact_charges) / len(intact_charges)
            intact_charges = [i - temp_correct for i in intact_charges]

        # Assign new charges to mol and then to atoms
        # NB: when the charges are assigned to atoms, degenerate atoms have an
        # averaged charge which is then reassigned to mol since mol is in atoms
        for index, atom in enumerate(mol):
            atom.q = intact_charges[index]
        assign_charges(mol, None, in_shell, None, max_bl)

        # New charges in mol
        intact_charges = [atom.q for atom in mol]

        # Damp the change in charges
        new_charges = [new * (1 - damping) + old * damping for new,
                       old in zip(intact_charges, old_charges)]
        # correct charges again (due to damping)
        if sum(new_charges) != 0.0:
            temp_correct = sum(new_charges) / len(new_charges)
            new_charges = [i - temp_correct for i in new_charges]
        # assign damped charges
        for index, atom in enumerate(mol):
            atom.q = new_charges[index]
        assign_charges(mol, None, in_shell, None, max_bl)

        # Calculate deviation between initial and new charges
        deviation = sum([abs(i - j)
                         for (i, j) in zip(intact_charges, old_charges)]) / len(mol)

        out_str = ("Iteration:", sc_loop, "Deviation:",
                   deviation, "Energy:", new_energy, "Charge self energy:", char_self, "Total - charge self:", new_energy - char_self)
        output_file.write(
            ("{:<6} {:<5} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f}\n".format(*out_str)))
        output_file.flush()

        return deviation

    def ewald_loop(in_atoms, in_mol, damping=0):
        """Do a single iteration of the Ewald + Gaussian loop

        Beware! This function changes the charge values of in_atoms and by
        extension in_mol. The point is to run this until convergence.

        Parameters:
        -----------
        in_atoms : list of Atom objects
            They should form a cell
        in_mol : list of Atom objects
            A list of atoms which are a part of in_atoms and form the quantum
            cluster
        damping : float
            Underrelaxation damping constant
        Returns:
        --------
        deviation : float
            Average deviation of charges from in_mol after an ewald calculation
            followed by a Gaussian calculation

        """
        global sc_loop
        sc_name = "sc_" + name
        sc_loop += 1

        # Initial charges in mol
        old_charges = [atom.q for atom in mol]
        # Ewald fitting
        run_ewald(sc_name, mol, in_atoms, vectors, in_nAt=nAt,
                  in_aN=aN, in_bN=bN, in_cN=cN, in_nChk=nChk)
        # read points output by Ewald
        sc_points = rf.read_points(sc_name + ".pts-cry")

        # Calculate new charges
        ef.write_gauss(sc_name + ".com", mol, sc_points,
                       os.path.join(here, sc_temp), proj_name=sc_name)
        subprocess.call("g09 " + sc_name + ".com", shell=True)

        intact_charges, new_energy, char_self, char_int = rf.read_g_char(
            sc_name + ".log", high_pop_method, debug=True)
        # Correct charges if they are not perfectly neutral
        if sum(intact_charges) != 0.0:
            # intact_charges[-1] -= sum(new_charges)
            temp_correct = sum(intact_charges) / len(intact_charges)
            intact_charges = [i - temp_correct for i in intact_charges]

        # Assign new charges to mol and then to atoms
        # NB: when the charges are assigned to atoms, degenerate atoms have an
        # averaged charge which is then reassigned to mol since mol is in atoms
        for index, atom in enumerate(mol):
            atom.q = intact_charges[index]
        assign_charges(mol, None, in_atoms, vectors, max_bl)

        # New charges in mol
        intact_charges = [atom.q for atom in mol]

        # Damp the change in charges
        new_charges = [new * (1 - damping) + old * damping for new,
                       old in zip(intact_charges, old_charges)]
        # correct charges again (due to damping)
        if sum(new_charges) != 0.0:
            temp_correct = sum(new_charges) / len(new_charges)
            new_charges = [i - temp_correct for i in new_charges]
        # assign damped charges
        for index, atom in enumerate(mol):
            atom.q = new_charges[index]
        assign_charges(mol, None, in_atoms, vectors, max_bl)

        # Calculate deviation between initial and new charges
        deviation = sum([abs(i - j)
                         for (i, j) in zip(intact_charges, old_charges)]) / len(mol)

        out_str = ("Iteration:", sc_loop, "Deviation:",
                   deviation, "Energy:", new_energy, "Charge self energy:", char_self, "Total - charge self:", new_energy - char_self)
        output_file.write(
            ("{:<6} {:<5} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f} {:<6} {:10.6f}\n".format(*out_str)))
        output_file.flush()

        return deviation

    output_file = open("prep.out", "w")

    # print start time
    start_time = datetime.now()
    output_file.write("STARTING TIME: " + str(start_time) + "\n")

    #-------------------------------------------------------------------------
    #-----------------------------READING INPUTS------------------------------
    #-------------------------------------------------------------------------

    # read config inputs
    inputs = pcf.parse_inputs("config")

    # name of the job goes here
    name = inputs["name"]

    # lattice vectors
    vectors = np.zeros((3, 3))
    # specified in config
    if all(len(i) == 3 for i in [inputs["a_vec"], inputs["b_vec"], inputs["c_vec"]]):
        a_vec = inputs["a_vec"]
        b_vec = inputs["b_vec"]
        c_vec = inputs["c_vec"]
        vectors[0] = a_vec
        vectors[1] = b_vec
        vectors[2] = c_vec
    else:  # from external file
        vectors = rf.read_vectors(inputs["vectors_file"])

    output_file.write("Vectors read in config:\n")
    output_file.write(str(vectors) + "\n")
    # name of the cell xyz file
    cell_file = inputs["cell_file"]

    # name of the xyz target shell file in case it is determined manually
    target_shell = inputs["target_shell"]

    # High level points specifications
    # name of the program for calculating charges
    high_pop_program = inputs["high_pop_program"]

    # name of the Gaussian log or cp2k file with population information
    high_pop_file = inputs["high_pop_file"]

    # kind of population in the Gaussian  or cp2k file: mulliken, esp or
    # hirshfeld
    high_pop_method = inputs["high_pop_method"]

    # Low level points specifications
    low_pop_program = inputs["low_pop_program"]
    low_pop_file = inputs["low_pop_file"]
    low_pop_method = inputs["low_pop_method"]

    # criterion to use for bond detection
    bonding = inputs["bonding"]

    # maximum bond length when defining a molecule
    bond_thresh = inputs["bond_thresh"]

    # label of an atom which will be part of the quantum cluster
    # warning: [0,N-1], not [1,N]
    atom_label = inputs["atom_label"]

    # the number of checkpoints in region 1
    nChk = inputs["nchk"]

    # the number of constrained charge atoms
    # i.e. atoms in regions 1 and 2
    nAt = inputs["nat"]

    # Ewald will multiply the unit cell in the direction
    # of the a, b or c vector 2N times (N positive and N negative)
    aN = inputs["an"]
    bN = inputs["bn"]
    cN = inputs["cn"]

    # the cluster will be of all molecules with atoms less than
    # clust_rad away from the centre of the central molecule
    clust_rad = inputs["clust_rad"]

    # how many times the input cluster needs to be repeated along each vector
    # positively and negatively to be able to contain the cluster to select.
    # the supercluster ends up being (1+2*traAN)*(1+2*traBN)*(1+2*traCN) times
    # bigger
    traAN = inputs["traan"]
    traBN = inputs["trabn"]
    traCN = inputs["tracn"]

    # use the self consistent version?
    self_consistent = inputs["self_consistent"]

    # Self Consistent Ewald Gaussian template
    sc_temp = inputs["sc_temp"]

    # Self Consistent Ewald deviation tolerance
    dev_tol = inputs["dev_tol"]

    # Ewald embedding
    ewald = inputs["ewald"]

    # Damping factor for underrelaxed sc loop. Use 0 for no Damping
    damping = inputs["damping"]

    # Print the cell with the complete model system at the origin?
    print_tweak = inputs["print_tweak"]

    # end config inputs

    #-------------------------------------------------------------------------
    #------------------------------END OF INPUTS------------------------------
    #-------------------------------------------------------------------------

    # read the input cell
    cell = rf.mol_from_file(cell_file)
    cell.vectors = vectors
    cell.bonding = bonding
    cell.thresh = bond_thresh
    atoms = rf.read_pos(cell_file)
    output_file.write("Read " + str(len(cell)) + " atoms in cell_file\n")
    output_file.flush()

    # High level charge assignment
    populate_cell(cell, high_pop_program, high_pop_file, high_pop_method)

    # get centered mod and corresponding cell
    model_system, mod_cell = cell.centered_mols(atom_label)
    # write useful xyz and new cell
    model_system.write_xyz("mol.init.xyz")
    if print_tweak:
        mod_cell.write_xyz("tweaked_cell.xyz")

    # Ewald section
    here = os.getcwd()

    if ewald:
        ew_path = os.path.join(here, 'ewald')
        if not os.path.exists(ew_path):
            os.makedirs(ew_path)
        os.chdir(ew_path)

        # Self Consistent EWALD
        if self_consistent:
            output_file.write("SELF CONSISTENT LOOP INITIATED\n")
            output_file.flush()
            sc_loop = 0
            while True:
                dev = ewald_loop(atoms, mol, damping)
                # check convergence
                if dev < dev_tol:
                    output_file.write("Tolerance reached: " +
                                      str(dev) + " < " + str(dev_tol) + "\n")
                    break

        # Final (or only) Ewald
        output_file.write("EWALD START\n")
        output_file.flush()
        run_ewald(name, mol, atoms, vectors, in_nAt=nAt,
                  in_aN=aN, in_bN=bN, in_cN=cN, in_nChk=nChk)
        output_file.write("EWALD END\n")
        output_file.flush()

        # read points output by Ewald
        points = rf.read_points(name + ".pts-cry")

    else:  # This means normal electrostatic embedding
        if target_shell:
            output_file.write("Reading the shell from: " + target_shell + "\n")
            output_file.flush()
            high_shell = rf.read_pos(target_shell)
            high_target_mol_char = rf.read_g_char(
                high_pop_file, high_pop_method)[0]
            # correct charges if they are not perfectly neutral
            if sum(high_target_mol_char) != 0.0:
                output_file.write("Charge correction: " +
                                  str(sum(high_target_mol_char)) + "\n")
                high_target_mol_char[-1] -= sum(high_target_mol_char)
            # assigns charges to a molecule
            high_target_pop_mol = rf.read_g_pos(high_pop_file)
            for index, atom in enumerate(high_target_pop_mol):
                atom.q = high_target_mol_char[index]

            # assign charges to the rest of the cell
            assign_charges(high_target_pop_mol, None, high_shell, None, max_bl)
        else:
            # get a cluster of atoms
            high_clust = Mol(atoms, vectors=vectors).make_cluster(clust_rad)
            # make a list of shell atoms
            high_shell = []
            for atom_i in high_clust:
                good_append = True
                for atom_j in mol:
                    if atom_i.very_close(atom_j):
                        good_append = False
                if good_append:
                    high_shell.append(atom_i)

        if self_consistent:
            output_file.write("SELF CONSISTENT LOOP INITIATED (no Ewald)\n")
            output_file.flush()
            sc_loop = 0
            while True:
                dev = ec_loop(high_shell, mol, damping)
                # check convergence
                if dev < dev_tol:
                    output_file.write("Tolerance reached: " +
                                      str(dev) + " < " + str(dev_tol) + "\n")
                    break

    os.chdir(here)
    # End of Ewald section

    # to manually input the cluster
    if target_shell:
        output_file.write("Reading the shell from: " + target_shell + "\n")
        shell = rf.read_pos(target_shell)
        target_mol_char = rf.read_g_char(low_pop_file, low_pop_method)[0]
        # correct charges if they are not perfectly neutral
        if sum(target_mol_char) != 0.0:
            output_file.write("Charge correction: " +
                              str(sum(target_mol_char)) + "\n")
            target_mol_char[-1] -= sum(target_mol_char)

        # assigns charges to a molecule
        target_pop_mol = rf.read_g_pos(low_pop_file)
        for index, atom in enumerate(target_pop_mol):
            atom.q = target_mol_char[index]

        # assign charges to the rest of the cell
        assign_charges(target_pop_mol, None, shell, None, max_bl)

    # to generate the cluster from radius
    else:
        # generate a shell of molecules with low level charges
        low_atoms = [copy(i) for i in atoms]
        populate_cell(low_atoms, low_pop_program, low_pop_file, low_pop_method)
        # get a cluster of atoms
        clust = Mol(low_atoms, vectors=vectors).make_cluster(clust_rad)
        # make a list of shell atoms
        shell = []
        for atom_i in clust:
            good_append = True
            for atom_j in mol:
                if atom_i.very_close(atom_j):
                    good_append = False
            if good_append:
                shell.append(atom_i)
        # write useful xyz
        ef.write_xyz("shell.xyz", shell)

    # Make inputs
    mh_path = os.path.join(here, 'mh')
    ml_path = os.path.join(here, 'ml')
    rl_path = os.path.join(here, 'rl')
    mg_path = os.path.join(here, 'mg')

    if not os.path.exists(mh_path):
        os.makedirs(mh_path)

    if not os.path.exists(ml_path):
        os.makedirs(ml_path)

    if not os.path.exists(rl_path):
        os.makedirs(rl_path)

    if not os.path.exists(mg_path):
        os.makedirs(mg_path)

    os.chdir(rl_path)
    ef.write_g_temp("rl.temp", shell, [], os.path.join(here, "rl.template"))
    os.chdir(ml_path)
    ef.write_g_temp("ml.temp", [], shell, os.path.join(here, "ml.template"))

    if ewald:
        os.chdir(mh_path)
        ef.write_g_temp("mh.temp", [], points,
                        os.path.join(here, "mh.template"))
        os.chdir(mg_path)
        ef.write_g_temp("mg.temp", [], points,
                        os.path.join(here, "mg.template"))  # only useful for CI
    else:
        # Make inputs
        os.chdir(mh_path)
        ef.write_g_temp("mh.temp", [], high_shell,
                        os.path.join(here, "mh.template"))
        os.chdir(mg_path)
        ef.write_g_temp("mg.temp", [], high_shell,
                        os.path.join(here, "mg.template"))  # only useful for CI
    os.chdir(here)
    end_time = datetime.now()
    output_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    output_file.write("ENDING TIME: " + str(end_time) + "\n")
    output_file.close()
