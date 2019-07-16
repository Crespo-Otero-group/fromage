#!/usr/bin/env python
"""This module is a tool to prepare template and input files for fromage.

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
import os
from datetime import datetime

from fromage.io import read_file as rf
from fromage.io import edit_file as ef
from fromage.io import parse_config_file as pcf
from fromage.scripts.fro_assign_charges import assign_charges
import fromage.utils.run_sequence as rs



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
        output_file = open(here+"/prep.out", "a")
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
            mol_char.bonding = in_mol.bonding
            mol_char.thresh = in_mol.thresh
            charges = [i.q for i in mol_char]
            # correct charges if they are not perfectly neutral
            if sum(charges) != 0.0:
                output_file.write("Charge correction: " +
                                  str(sum(charges)) + "\n")
                mol_char[-1].q -= sum(charges)

            # assign charges to the rest of the cell
            in_mol.populate(mol_char)

        output_file.close()
        return

    here = os.getcwd()
    output_file = open(here+"/prep.out", "w")

    # print start time
    start_time = datetime.now()
    output_file.write("STARTING TIME: " + str(start_time) + "\n")

    # read config inputs
    inputs = pcf.parse_inputs("config")

    # read the input cell
    cell = rf.mol_from_file(inputs["cell_file"])
    cell.vectors = inputs["vectors"]
    cell.bonding = inputs["bonding"]
    cell.thresh = inputs["bond_thresh"]
    cell = cell.confined()

    output_file.write("Read " + str(len(cell)) + " atoms in cell_file\n")
    output_file.close()

    # High level charge assignment
    populate_cell(cell, inputs["high_pop_program"], inputs["high_pop_file"], inputs["high_pop_method"])
    region_1, cell = cell.centered_mols(inputs["atom_label"])

    # write useful xyz and new cell
    region_1.write_xyz("mol.init.xyz")
    if inputs["print_tweak"]:
        ef.write_xyz("tweaked_cell.xyz", cell)

    run_sequence = rs.RunSeq(region_1, cell, inputs)

    region_2, high_points = run_sequence.run()

    region_2.write_xyz("shell.xyz")

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
    ef.write_g_temp("rl.temp", region_2, [], os.path.join(here, "rl.template"))
    os.chdir(ml_path)
    ef.write_g_temp("ml.temp", [], region_2, os.path.join(here, "ml.template"))

    os.chdir(mh_path)
    ef.write_g_temp("mh.temp", [], high_points,
                    os.path.join(here, "mh.template"))
    os.chdir(mg_path)
    ef.write_g_temp("mg.temp", [], high_points,
                    os.path.join(here, "mg.template"))
    os.chdir(here)
    end_time = datetime.now()

    output_file = open(here+"/prep.out", "a")
    output_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
    output_file.write("ENDING TIME: " + str(end_time) + "\n")
    output_file.close()
