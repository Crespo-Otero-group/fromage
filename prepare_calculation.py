#!/usr/bin/env python
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
from assign_charges import assign_charges
from atom import Atom
from datetime import datetime

    def run_ewald(in_name,in_mol,in_atoms,in_vectors,in_nAt=500,in_aN=2,in_bN=2,in_cN=2,in_nChk=1000):
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
        ef.write_uc(in_name + ".uc", in_vectors, in_aN, in_bN, in_cN, in_atoms)
        ef.write_qc(in_name + ".qc", in_mol)
        ef.write_ew_in(in_name, "ewald.in." + in_name, in_nChk, in_nAt)
        ef.write_seed()
        # run Ewald
        subprocess.call("./Ewald < ewald.in." + in_name, shell=True)

if __name__ == '__main__':
    output_file = open("prep.out", "w")

    # print start time
    start_time = datetime.now()
    output_file.write("STARTING TIME: " + str(start_time) + "\n")

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

    output_file.write("Vectors read in config:\n")
    output_file.write(str(vectors) + "\n")
    # name of the cell xyz file
    if "cell_file" in inputs:
        cell_file = inputs["cell_file"]
    else:
        cell_file = name + "xyz"

    # name of the cp2k file with population information
    if "cp2k_file" in inputs:
        cp2k_file = inputs["cp2k_file"]
    else:
        cp2k_file = ""

    # name of the Gaussian log file with population information
    if "mol_pop_file" in inputs:
        mol_pop_file = inputs["mol_pop_file"]
    else:
        mol_pop_file = ""

    # kind of population in the Gaussian file 0:Mulliken 1:ESP
    if "mol_pop_kind" in inputs:
        mol_pop_kind = int(inputs["mol_pop_kind"])
    else:
        mol_pop_kind = 0

    # maximum bond length when defining a molecule
    if "max_bl" in inputs:
        max_bl = float(inputs["max_bl"])
    else:
        max_bl = 1.7

    # label of an atom which will be part of the quantum cluster
    # warning: [0,N-1], not [1,N]
    if "label_atom" in inputs:
        if type(inputs["label_atom"]) == str:
            label_atom = [int(inputs["label_atom"])]
        else:
            label_atom = [int(i) for i in inputs["label_atom"]]
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
    if "cp2k_pop_method" in inputs:
        cp2k_pop_method = int(inputs["cp2k_pop_method"])
    else:
        cp2k_pop_method = 0

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

    # Self Consistent Ewald Gaussian template
    if "sc_temp" in inputs:
        sc_temp = inputs["sc_temp"]
    else:
        sc_temp = ""

    # Self Consistent Ewald atom kind 0:Mulliken 1:ESP
    if "sc_kind" in inputs:
        sc_kind = int(inputs["sc_kind"])
    else:
        sc_kind = 0

    # Self Consistent Ewald deviation tolerance
    if "dev_tol" in inputs:
        dev_tol = float(inputs["dev_tol"])
    else:
        dev_tol = 0.001

    # Ewald embedding, use 0 for false
    if "ewe" in inputs:
        ewe = bool(inputs["ewe"])
    else:
        ewe = True

    # Cluster self consistent template for mol
    if "csc_temp_h" in inputs:
        csc_temp_h = inputs["csc_temp_h"]
    else:
        csc_temp_h = ""

    # Cluster self consistent template for shell
    if "csc_temp_l" in inputs:
        csc_temp_l = inputs["csc_temp_l"]
    else:
        csc_temp_l = ""

    # end config inputs


    # read the input atoms
    atoms = rf.read_xyz(cell_file)[-1]
    output_file.write("Read " + str(len(atoms)) + " atoms in cell_file\n")

    # the molecule of interest and the atoms which now contain
    # the full, unchopped molecule
    # NB: all objects in mol are also referenced inside atoms
    mol, atoms = ha.complete_mol(max_bl, atoms, label_atom, vectors)

    # read charges
    if cp2k_file:
        charges = rf.read_cp2k(cp2k_file, cp2k_pop_method)[0]
        output_file.write("Read " + str(len(atoms)) + " charges in cp2k_file\n")
        # in case there are more charges than atoms for some reason
        charges = charges[:len(atoms)]
        # correct charges if they are not perfectly neutral
        if sum(charges) != 0.0:
            output_file.write("Charge correction: " + str(sum(charges)) + "\n")
            charges[-1] -= sum(charges)

        # assigns charges to atoms
        for index, atom in enumerate(atoms):
            atom.q = charges[index]
    elif mol_pop_file:
        mol_char = rf.read_g_char(mol_pop_file, mol_pop_kind)[0]
        # correct charges if they are not perfectly neutral
        if sum(mol_char) != 0.0:
            output_file.write("Charge correction: " + str(sum(mol_char)) + "\n")
            mol_char[-1] -= sum(mol_char)

        # assigns charges to a molecule
        pop_mol = rf.read_g_pos(mol_pop_file)
        for index, atom in enumerate(pop_mol):
            atom.q = mol_char[index]

        # assign charges to the rest of the cell
        assign_charges(pop_mol, None, atoms, vectors, max_bl)
    else:
        output_file.write(
            "Please input a population file with cp2k_file or mol_pop_file\n")

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

    # Self Consistent EWALD
    #def loop_ewald(in_name,in_mol,in_vectors,in_atoms,in_aN,in_bN,in_cN,in_nChk,in_nAt,in_sc_kind,in_max_bl,in_dev_tol):
    if sc_temp:
        output_file.write("SELF CONSISTENT LOOP INITIATED\n")
        sc_name = "sc_" + name
        sc_loop = 0
        while True:
            sc_loop += 1
            old_charges = [atom.q for atom in mol]
            # Ewald fitting
            run_ewald(sc_name,mol,atoms,vectors,nat,aN,bN,cN,nChk)
            # read points output by Ewald
            sc_points = rf.read_points(sc_name + ".pts-tb")

            ef.write_gauss(sc_name, sc_name + ".com", mol, sc_points, sc_temp)
            # calculate new charges
            subprocess.call("g09 " + sc_name + ".com", shell=True)
            new_charges = rf.read_g_char(sc_name + ".log", sc_kind)[0]

            # correct charges if they are not perfectly neutral
            if sum(new_charges) != 0.0:
                new_charges[-1] -= sum(new_charges)

            # assign new charges
            for index, atom in enumerate(mol):
                atom.q = new_charges[index]
            assign_charges(mol, None, atoms, vectors, max_bl)
            # spread evenly
            new_charges = [atom.q for atom in mol]

            # assign new charges
            for index, atom in enumerate(mol):
                atom.q = new_charges[index]
            assign_charges(mol, None, atoms, vectors, max_bl)

            # check convergence
            deviation = sum([abs(i - j)
                             for (i, j) in zip(new_charges, old_charges)]) / len(mol)
            output_file.write("Iteration: " + str(sc_loop) +
                              "  Deviation: " + str(deviation) + "\n")
            if deviation < dev_tol:
                output_file.write("Tolerance reached: " +
                                  str(deviation) + " < " + str(dev_tol) + "\n")
                break
            output_file.flush()
            # assign new charges
            for index, atom in enumerate(mol):
                atom.q = new_charges[index]
            assign_charges(mol, None, atoms, vectors, max_bl)

    #Self consistent between cluster and mol, no Ewald
    elif csc_temp_h:
        csc_name_h = "csc_" + name + "_h"
        csc_name_l = "csc_" + name + "_l"
        csc_loop = 0
        while True:
            csc_loop += 1
            old_charges_h = [atom.q for atom in mol]
            old_charges_l = [atom.q for atom in shell]

            # calculate the middle molecule embedded in points
            ef.write_gauss(csc_name_h, csc_name_h + ".com", mol, shell, csc_temp_h)
            subprocess.call("g09 " + csc_name_h + ".com", sell=True)
            new_charges_h = rf.read_g_char(csc_name_h + ".log", sc_kind)[0]

            # correct for neutrality
            if sum(new_charges) != 0.0:
                new_charges[-1] -= sum(new_charges)

            # calculate the surrounding shell
            ef.write_gauss(csc_name_l, csc_name_l + ".com", shell, mol, csc_temp_l)
            subprocess.call("g09 " + csc_name_l + ".com", sell=True)
            new_charges_l = rf.read_g_char(csc_name_l + ".log", sc_kind)[0]

            # assign new charges
            for index, atom in enumerate(mol):
                atom.q = new_charges_h[index]

            for index, atom in enumerate(shell):
                atom.q = new_charges_l[index]

            # check convergence
            deviation_h = sum([abs(i - j)
                               for (i, j) in zip(new_charges_h, old_charges_h)]) / len(mol)
            deviation_l = sum([abs(i - j)
                               for (i, j) in zip(new_charges_l, old_charges_l)]) / len(shell)
            output_file.write("Iteration: " + str(csc_loop) +
                              "\nDeviation_h: " + str(deviation_h) +
                              "\nDeviation_l: " + str(deviation_l))
            if deviation_h < dev_tol and deviation_l < dev_tol:
                output_file.write("Tolerance reached:" + str(deviation_h) +
                                  " and " + str(deviation_l) + " < " + dev_tol)
                break
            output_file.flush()

    # Final (or only) Ewald
    if ewe:
        # write inputs
        ef.write_uc(name + ".uc", vectors, aN, bN, cN, atoms)
        ef.write_qc(name + ".qc", mol)
        ef.write_ew_in(name, "ewald.in." + name, nChk, nAt)
        ef.write_seed()
        # run Ewald
        subprocess.call("./Ewald < ewald.in." + name, shell=True)
        # read points output by Ewald
        points = rf.read_points(name + ".pts-tb")

    if ewe:
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
        ef.write_g_temp("mh", "mh.temp", [], shell, "mh.template")
        ef.write_g_temp("mg", "mg.temp", [], shell,
                        "mg.template")  # only useful for CI
    end_time = datetime.now()
    output_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    output_file.write("ENDING TIME: " + str(end_time) + "\n")
    output_file.close()
