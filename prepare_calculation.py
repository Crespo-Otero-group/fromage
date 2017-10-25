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
import read_file as rf
import edit_file as ef
import handle_atoms as ha
import parse_config_file as pcf
from assign_charges import assign_charges
from atom import Atom
from datetime import datetime
from copy import copy


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
    def populate_cell(in_atoms, program, pop_file, method):
        """
        Assign charge to the atoms from a unit cell

        Don't import this, it depends on variables in the main of this module:
        output_file, vectors and max_bl
        If you want to assign charges go straign to assign_charges.py and use
        those utilities.

        Parameters
        ----------
        in_atoms : list of Atom objects
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
            output_file.write("Read " + str(len(in_atoms)) +
                              " charges in cp2k_file\n")
            # in case there are more charges than atoms
            charges = charges[:len(in_atoms)]
            # correct charges if they are not perfectly neutral
            if sum(charges) != 0.0:
                output_file.write("Charge correction: " +
                                  str(sum(charges)) + "\n")
                charges[-1] -= sum(charges)

        if program.lower() == "gaussian":
            mol_char = rf.read_g_char(pop_file, method)[0]
            # correct charges if they are not perfectly neutral
            if sum(mol_char) != 0.0:
                output_file.write("Charge correction: " +
                                  str(sum(mol_char)) + "\n")
                mol_char[-1] -= sum(mol_char)

            # assigns charges to a molecule
            pop_mol = rf.read_g_pos(pop_file)
            for index, atom in enumerate(pop_mol):
                atom.q = mol_char[index]

            # assign charges to the rest of the cell
            assign_charges(pop_mol, None, in_atoms, vectors, max_bl)
        return

    output_file = open("prep.out", "w")

    # print start time
    start_time = datetime.now()
    output_file.write("STARTING TIME: " + str(start_time) + "\n")

    def ewald_loop(in_atoms, in_mol):
        """A single iteration of the Ewald + Gaussian loop

        Beware! This function changes the charge values of in_atoms and by
        extension in_mol. The point is to run this until convergence.

        Parameters:
        -----------
        in_atoms : list of Atom objects
            They should form a cell
        in_mol: list of Atom objects
            A list of atoms which are a part of in_atoms and form the quantum
            cluster
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
        run_ewald(sc_name, mol, atoms, vectors, in_nAt=nAt,
                  in_aN=aN, in_bN=bN, in_cN=cN, in_nChk=nChk)
        # read points output by Ewald
        sc_points = rf.read_points(sc_name + ".pts-tb")

        # Calculate new charges
        ef.write_gauss(sc_name, sc_name + ".com", mol, sc_points, sc_temp)
        subprocess.call("g09 " + sc_name + ".com", shell=True)
        new_charges = rf.read_g_char(sc_name + ".log", high_pop_method)[0]

        # Correct charges if they are not perfectly neutral
        if sum(new_charges) != 0.0:
            new_charges[-1] -= sum(new_charges)

        # Assign new charges to mol and then to atoms
        # NB: when the charges are assigned to atoms, degenerate atoms have an
        # averaged charge which is then reassigned to mol since mol is in atoms
        for index, atom in enumerate(mol):
            atom.q = new_charges[index]
        assign_charges(mol, None, atoms, vectors, max_bl)

        # New charges in mol
        new_charges = [atom.q for atom in mol]

        # Calculate deviation between initial and new charges
        deviation = sum([abs(i - j)
                         for (i, j) in zip(new_charges, old_charges)]) / len(mol)
        output_file.write("Iteration: " + str(sc_loop) +
                          "  Deviation: " + str(deviation) + "\n")
        output_file.flush()

        return deviation

    #-------------------------------------------------------------------------
    #-----------------------------READING INPUTS------------------------------
    #-------------------------------------------------------------------------

    # read config inputs
    inputs = pcf.parse_inputs("config")

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

    output_file.write("Vectors read in config:\n")
    output_file.write(str(vectors) + "\n")
    # name of the cell xyz file
    cell_file = inputs["cell_file"]

    # name of the xyz target shell file in case it is determined manually
    target_clust = inputs["target_clust"]

    # High level points specifications
    # name of the program for calculating charges
    high_pop_program = inputs["high_pop_program"]

    # name of the Gaussian log or cp2k file with population information
    high_pop_file = inputs["high_pop_file"]

    # kind of population in the Gaussian  or cp2k file 0:Mulliken 1:ESP 2:
    # Hirshfeld
    high_pop_method = inputs["high_pop_method"]

    # Low level points specifications
    low_pop_program = inputs["low_pop_program"]
    low_pop_file = inputs["low_pop_file"]
    low_pop_method = inputs["low_pop_method"]

    # maximum bond length when defining a molecule
    max_bl = float(inputs["max_bl"])

    # label of an atom which will be part of the quantum cluster
    # warning: [0,N-1], not [1,N]
    atom_label = inputs["atom_label"]

    # the number of checkpoints in region 1
    nChk = int(inputs["nchk"])

    # the number of constrained charge atoms
    # i.e. atoms in regions 1 and 2
    nAt = int(inputs["nat"])

    # Ewald will multiply the unit cell in the direction
    # of the a, b or c vector 2N times (N positive and N negative)
    aN = int(inputs["an"])
    bN = int(inputs["bn"])
    cN = int(inputs["cn"])

    # the cluster will be of all molecules with atoms less than
    # clust_rad away from the centre of the central molecule
    clust_rad = float(inputs["clust_rad"])

    # how many times the input cluster needs to be repeated along each vector
    # positively and negatively to be able to contain the cluster to select.
    # the supercluster ends up being (1+2*traAN)*(1+2*traBN)*(1+2*traCN) times
    # bigger
    traAN = int(inputs["traan"])
    traBN = int(inputs["trabn"])
    traCN = int(inputs["tracn"])

    # use the self consistent version?
    self_consistent = inputs["self_consistent"]

    # Self Consistent Ewald Gaussian template
    sc_temp = inputs["sc_temp"]

    # Self Consistent Ewald deviation tolerance
    dev_tol = float(inputs["dev_tol"])

    # Ewald embedding, use 0 for false
    ewald = bool(inputs["ewald"])

    # end config inputs

    #-------------------------------------------------------------------------
    #------------------------------END OF INPUTS------------------------------
    #-------------------------------------------------------------------------

    # read the input atoms
    atoms = rf.read_pos(cell_file)
    output_file.write("Read " + str(len(atoms)) + " atoms in cell_file\n")

    # the molecule of interest and the atoms which now contain
    # the full, unchopped molecule
    # NB: all objects in mol are also referenced inside atoms
    mol, atoms = ha.complete_mol(max_bl, atoms, atom_label, vectors)

    # High level charge assignment
    populate_cell(atoms, high_pop_program, high_pop_file, high_pop_method)

    # find the centroid of the molecule
    c_x, c_y, c_z = ha.find_centroid(mol)
    # translate the molecule and atoms to the centroid
    for atom in atoms:
        atom.translate(-c_x, -c_y, -c_z)

    # write useful xyz
    ef.write_xyz("mol.init.xyz", mol)

    # make a very big cell
    mega = ha.make_mega_cell(atoms, traAN, traBN, traCN, vectors)
    # get a cluster of atoms
    clust = ha.make_cluster(mega, clust_rad, max_bl)

    # Self Consistent EWALD
    if self_consistent:
        output_file.write("SELF CONSISTENT LOOP INITIATED\n")
        sc_loop = 0
        while True:
            dev = ewald_loop(atoms, mol)
            # check convergence
            if dev < dev_tol:
                output_file.write("Tolerance reached: " +
                                  str(dev) + " < " + str(dev_tol) + "\n")
                break

    # Final (or only) Ewald
    if ewald:
        run_ewald(name, mol, atoms, vectors, in_nAt=nAt,
                  in_aN=aN, in_bN=bN, in_cN=cN, in_nChk=nChk)
        # read points output by Ewald
        points = rf.read_points(name + ".pts-tb")

    else:  # This means normal electrostatic embedding
        if target_clust:
                out_file.write("Reading the shell from: "+target_clust+"\n")
                high_shell = rf.read_pos(target_clust)
                high_target_mol_char = rf.read_g_char(high_pop_file, high_pop_method)[0]
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
                assign_charges(target_pop_mol, None, shell, None, max_bl)
        else:
            # make a very big cell
            high_mega = ha.make_mega_cell(atoms, traAN, traBN, traCN, vectors)
            # get a cluster of atoms
            high_clust = ha.make_cluster(high_mega, clust_rad, max_bl)
            # make a list of shell atoms
            high_shell = []
            for atom_i in high_clust:
                good_append = True
                for atom_j in mol:
                    if atom_i.very_close(atom_j):
                        good_append = False
                if good_append:
                    high_shell.append(atom_i)

    # to manually input the cluster
    if target_clust:
        out_file.write("Reading the shell from: "+target_clust+"\n")
        shell = rf.read_pos(target_clust)
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
        # make a very big cell
        mega = ha.make_mega_cell(low_atoms, traAN, traBN, traCN, vectors)
        # get a cluster of atoms
        clust = ha.make_cluster(mega, clust_rad, max_bl)
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
        ef.write_xyz("clust.xyz", clust)
        ef.write_xyz("shell.xyz", shell)


    if ewald:
        # Make inputs
        ef.write_g_temp("rl", "rl.temp", shell, [], "rl.template")
        ef.write_g_temp("ml", "ml.temp", [], shell, "ml.template")
        ef.write_g_temp("mh", "mh.temp", [], points, "mh.template")
        ef.write_g_temp("mg", "mg.temp", [], points,
                        "mg.template")  # only useful for CI
    else:
        # Make inputs
        ef.write_g_temp("rl", "rl.temp", shell, [], "rl.template")
        ef.write_g_temp("ml", "ml.temp", [], shell, "ml.template")
        ef.write_g_temp("mh", "mh.temp", [], high_shell, "mh.template")
        ef.write_g_temp("mg", "mg.temp", [], high_shell,
                        "mg.template")  # only useful for CI
    end_time = datetime.now()
    output_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    output_file.write("ENDING TIME: " + str(end_time) + "\n")
    output_file.close()
