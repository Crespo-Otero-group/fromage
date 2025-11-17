#!/usr/bin/env python
"""
Script to generate point-charge embedded Gaussian16 calculations for diabatisation with Overdia.

Currently only for agg but can be extedended to aggregates. Overdia must use TDDFT or TDA.

Michael Ingham 16/04/24
"""
import os
from datetime import datetime
import numpy as np 

from fromage.io import read_file as rf
from fromage.io import edit_file as ef
from fromage.io import parse_config_file as pcf
import fromage.utils.run_sequence as rs
import fromage.scripts.fro_assign_charges as fc
from fromage.utils.atom import Atom

from fromage.utils.mol import Mol
import argparse
import subprocess

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
    outfile = open(os.getcwd() + "/fro_overdia.out", "a")
    if program.lower() == "cp2k":
        charges = rf.read_cp2k(pop_file, method)[0]
        outfile.write("Read " + str(len(in_mol)) + " charges in cp2k_file\n")
        # in case there are more charges than atoms
        charges = charges[: len(in_mol)]
        # correct charges if they are not perfectly neutral
        if sum(charges) != 0.0:
            outfile.write("Charge correction: " + str(sum(charges)) + "\n")
            charges[-1] -= sum(charges)

    if program.lower() == "gaussian":
        mol_char = rf.mol_from_gauss(pop_file, pop=method)
        print(mol_char)
        mol_char.bonding = in_mol.bonding
        mol_char.thresh = in_mol.thresh
        charges = [i.q for i in mol_char]
        print(mol_char.bonding, mol_char.thresh)
        # correct charges if they are not perfectly neutral
        if sum(charges) != 0.0:
            outfile.write("Charge correction: " + str(sum(charges)) + "\n")
            mol_char[-1].q -= sum(charges)

        # assign charges to the rest of the cell
        print(in_mol.bonding, in_mol.thresh)
        in_mol.populate(mol_char)

    outfile.close()
    return


def getPrograms():
    """
    Read requested programs for high_level and low_level calculations
        from fromage.in file to determine which template files to write

    Returns
    -------
    high_level_write : Function object from io.edit_file.py (ef) that
        writes the correct template file for the high level method
    low_level_write : Function object from io.edit_file.py (ef) that
        writes the correct template file for the low level method

    """
    def_inputs = {"high_level": "gaussian", "low_level": "gaussian"}

    inputs = def_inputs.copy()

    if os.path.isfile("fromage.in"):
        new_inputs = rf.read_config("fromage.in")
        inputs.update(new_inputs)

    writer_list = []
    for prog in [inputs["high_level"], inputs["low_level"]]:
        if prog == "xtb":
            writer_list.append(ef.write_xtb_temp)
        elif prog == "fomo-ci" or prog == "mopac":
            writer_list.append(ef.write_tinker_temp)
        else:
            writer_list.append(ef.write_g_temp)
    return writer_list[0], writer_list[1]


def neutralise_cluster(clust, outfile):
    """Remove any net charge from cluster"""
    net_charge = sum([atom.q for atom in clust])

    #  neutralise 
    corr = net_charge / len(clust)
    for atom in clust:
        atom.q -= corr
    net_charge = sum([atom.q for atom in clust])
    outfile.write(f"\nNet charge: {net_charge:.2e} (Correction: {corr:.2e}) ")
   
    return clust


def set_gauss(filep, nstates, type, nprocs):
    """
    run sed command to add checkpoint string and number of excited states for each. Probably can do this in a better way
    """
    cmd = f"sed  -i 's/xxxnstatesxxx/{str(nstates)}/g' {filep}"
    subprocess.run(cmd, shell=True)

    cmd = f"sed  -i 's/xxxchkxxx/{type}/g' {filep}"
    subprocess.run(cmd, shell=True)

    cmd = f"sed  -i 's/xxxnprocsxxx/{str(nprocs)}/g' {filep}"
    subprocess.run(cmd, shell=True)
    return


def vis_pce(mol, charges, name, outfile=None):
    """
    change point atom.elem to hydrogen and add to molecule, then write xyz.
    """
    pc_path = "vis/"
    if not os.path.exists(pc_path):
        os.mkdir(pc_path)

    vis = mol.copy()
    for char in charges:
        char_vis = char.copy()
        char_vis.elem = "H"
        vis.append(char_vis)

    vis.remove_duplicates()
    if outfile:
        outfile.write("Number of point charges: ", len(vis))
    vis.write_xyz(name)
    return

def main(
    n_agg_states,
    n_mono_states,
    nprocs,
    vis_charges,
    run_type,
    parallel_job,
    custom_region,
    reuse_charges
):
    """
    Run RunSequnce object to get point charges, then create three files: mono1.com, mono2.com and agg

    So far, Ewald is performed for just the agg, then the local charges modified. Indices is used to specify the fragments for some other scripts I've written
    """
    here = os.getcwd()

    if os.path.exists("fro_overdia.out"):
        os.remove("fro_overdia.out")

    outfile = open(here + "/fro_overdia.out", "w")

    # print start time
    start_time = datetime.now()
    outfile.write("STARTING TIME: " + str(start_time) + "\n")

    # read config inputs
    inputs = pcf.parse_inputs("config")

    # read the input cell
    cell = rf.mol_from_file(inputs["cell_file"])
    cell.vectors = inputs["vectors"]
    cell.bonding = inputs["bonding"]
    cell.thresh = inputs["bond_thresh"]
    cell = cell.confined()

    # write some nice output
    outfile.write("-------------------------\n     Input information\n-------------------------")
    outfile.write(f"\nAtoms in unit cell: {len(cell)}")
    outfile.write(f"\nBonding: {cell.thresh}-{cell.bonding}")

    # get indices from config
    indices = inputs["atom_label"]
    outfile.write("\nAtom labels: " + " ".join(map(str, np.array(indices)+1)))
    if len(indices) < 2:
        raise ValueError(
            "At least two molecules must be included in model for fragment diabatisation"
        )
    n_monomers = len(indices)

    # High level charge assignment to cell
    populate_cell(
        cell,
        inputs["high_pop_program"],
        inputs["high_pop_file"],
        inputs["high_pop_method"],
    )
    outfile.write("\nPopulation file: " + inputs["high_pop_file"])
    region_1, cell = cell.centered_mols(indices)  ## use fragments specified by argparse
    region_1_pc = region_1.copy()

    region_1.write_xyz("mol.init.xyz")
    outfile.write("\nWriting mol.init.xyz")
    if inputs["print_tweak"]:
        ef.write_xyz("tweaked_cell.xyz", cell)


    # check len of region_1 is the same as expected
    n_fragments = len(region_1.segregate())
    if not n_fragments == n_monomers:
        raise ValueError("Different number of monomers than specified in config")
    else:
        outfile.write(f"\nFragments detected: {n_fragments}")

    # read in charge distribution from previous run
    if reuse_charges:
        outfile.write("\n----------------------\n     Restart\n----------------------\n")
        region_2 = rf.mol_from_file("shell.xyz")
        outfile.write("\nReading in shell.xyz")
        high_points = Mol([])

        outfile.write("\nReading in high_charges.pc")
        with open("high_charges.pc", "r") as pc_file:
            next(pc_file)
            for pc_line in pc_file:
                pc_line =  pc_line.split()
                pc = Atom("H", float(pc_line[0]), float(pc_line[1]), float(pc_line[2]), float(pc_line[3]))
                # print(pc)
                high_points.append(pc)
        
        outfile.write("\nReading in low_charges.pc")
        with open("low_charges.pc", "r") as pc_file:
            next(pc_file)
            for i, pc in enumerate(pc_file):
                region_2[i].q = float(pc.split()[-1])

    else:
        # generate shell region and PCE
        outfile.write("\n-------------------------\n     Calling RunSeq\n-------------------------\n")
        outfile.close()
        run_sequence = rs.RunSeq(region_1, cell, inputs, out_file_name="fro_overdia.out")
        region_2, high_points = run_sequence.run()
        region_2.write_xyz("shell.xyz")
        
        #reopen output
        outfile = open(here + "/fro_overdia.out", "a")

        ## store point charges in files for reruns
        outfile.write("Saving charges to high_charges.pc")
        with open("high_charges.pc", "w") as pc_file:
            pc_file.write("{}\n".format(len(high_points)))
            for pc in high_points:
                pc_file.write("{} {} {} {}\n".format(
                pc.x, pc.y, pc.z, pc.q
                ))
        
        outfile.write("Saving charges to low_charges.pc")
        with open("low_charges.pc", "w") as pc_file:
            pc_file.write("{}\n".format(len(region_2)))
            for pc in region_2:
                pc_file.write("{} {} {} {}\n".format(
                pc.x, pc.y, pc.z, pc.q
                ))


        outfile = open(here + "/fro_overdia.out", "a")
        outfile.write("Writing shell.xyz")


    outfile.write("\n-------------------------------------\n     Generating FrD(EE) input\n-------------------------------------")
    ## overwrite generated file read-in custom QM region (e.g. defect or optimised)
    if custom_region != None:
        outfile.write("\nModel region updated from model.init.xyz")
        region_1 = rf.mol_from_file("model.init.xyz")
        fc.assign_charges(region_1_pc, region_1)

    # neutralise shell and write agg file
    outfile.write(f"\nEnsuring fromage point charges are neutral")
    high_points = neutralise_cluster(high_points, outfile)

    # divide agg into respective monomers and reconsitute; consistency of atom indexing between (mono1 + mono2) and agg
    # is required for Overdia to run
    agg_as_mols = region_1.segregate(diff_mols=False)
    outfile.write(f"\nReordering aggregate atoms by monomer")

    ## iterate to get list of fragments and the aggregate object. NB: the atoms in fragments and monomers *must* be in same order
    monomers, mono_envs, mono_paths, agg = [], [], [], Mol([])

    # iterate through number of fragments
    for i in range(n_monomers):

        # get fragments
        monomers.append(agg_as_mols[i])

        # reorder agg to match monomers
        agg.extend(agg_as_mols[i])

        # get list of environment cops
        outfile.write(f"\nGetting monomer {i+1} charges")
        mono_envs.append(neutralise_cluster(high_points.copy(), outfile))

        # get file path
        mono_paths.append(f"mono{i+1}.com")

    # Make aggregate gaussian16
    outfile.write(f"\n\nAggregate")
    outfile.write(f"\nWriting G16 file 'agg.com'")
    agg_p = "agg.com"
    ef.write_gauss(agg_p, agg, high_points, "gauss.temp")
    set_gauss(agg_p, nstates=n_agg_states, type="agg", nprocs=nprocs)

    # visualise
    if vis_charges:
        outfile.write(f"\nGenerating visualisation 'vis/agg.xyz'")
        vis_pce(agg, high_points, "vis/agg.xyz")

    # make the gaussian16 input
    for i, (mono, mono_env, path) in enumerate(zip(monomers, mono_envs, mono_paths)):
        outfile.write(f"\n\nMonmer {i+1}")
        ## add missing point charges from dimer calculation
        for j in range(len(agg_as_mols)):
            if j != i:
                outfile.write(f"\nAdding charges to monomer environment")
                mono_env.extend(agg_as_mols[j])

        # write gaussian input
        outfile.write(f"\nWriting G16 file 'mono{i}.com'")
        ef.write_gauss(path, mono, mono_env, "gauss.temp")
        set_gauss(path, nstates=n_mono_states, type=f"mono{i+1}", nprocs=nprocs)

        # visualise
        if vis_charges:
            outfile.write(f"\nGenerating visualisation 'vis/mono{i+1}.xyz'")
            vis_pce(mono, mono_env, f"vis/mono{i+1}.xyz")

    end_time = datetime.now()
    outfile.write("\n\nELAPSED TIME: " + str(end_time - start_time) + "\n")
    outfile.write("ENDING TIME: " + str(end_time) + "\n")
    outfile.close()
    outfile.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A script for generating G16 input files (i.e. mono1.com, mono2.com, and agg.com) for Overdia diabatisation calculation. The customizable G16 template file (gauss.temp) must be provided."
    )

    # Define optional arguments for variable settings
    parser.add_argument(
        "--N",
        type=int,
        default=60,
        help="N aggregate excited states",
    )
    parser.add_argument(
        "--M",
        type=int,
        default=10,
        help="M monomer excited states",
    )
    parser.add_argument(
        "--nprocs", type=int, default=40, help="nproc processors in input file"
    )
    parser.add_argument(
        "--vis_charges", type=bool, default=True, help="Visualize the electrostatic embedding as .xyz files written to vis/"
    )
    parser.add_argument(
        "--run_type",
        type=str,
        default="sp",
        help="sp: generate files for single-point FrD(EE), freq: generate input for normal mode calculations",
    )
    parser.add_argument(
        "--parallel_job", type=bool, default=True, help="For FrD-LVC(EE): Generate displacements in parallel"
    )
    parser.add_argument(
        "--custom_region", type=str, default=None, help="Manually specify mol.init.xyz from model.init.xyz file. Useful for investigating different parts of the PES or defects."
    )
    parser.add_argument('--reuse_charges', type=str, default=False, help='Read charges.pc file, Requires mol.init.xyz, shell.xyz as well')


    args = parser.parse_args()
    main(
        args.N,
        args.M,
        args.nprocs,
        args.vis_charges,
        args.run_type,
        args.parallel_job,
        args.custom_region,
        args.reuse_charges
    )
