#!/usr/bin/env python
"""
Script to generate point-charge embedded Gaussian16 calculations for diabatisation with Overdia.

Currently only for agg but can be extedended to aggregates. Overdia must use TDDFT or TDA.

Michael Ingham 16/04/24
"""
import os
from datetime import datetime

from fromage.io import read_file as rf
from fromage.io import edit_file as ef
from fromage.io import parse_config_file as pcf
import fromage.utils.run_sequence as rs
import fromage.scripts.fro_assign_charges as fc

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
    outfile= open(os.getcwd() + "/frooverdia.out", "a")
    if program.lower() == "cp2k":
        charges = rf.read_cp2k(pop_file, method)[0]
        outfile.write("Read " + str(len(in_mol)) +
                            " charges in cp2k_file\n")
        # in case there are more charges than atoms
        charges = charges[:len(in_mol)]
        # correct charges if they are not perfectly neutral
        if sum(charges) != 0.0:
            outfile.write("Charge correction: " +
                                str(sum(charges)) + "\n")
            charges[-1] -= sum(charges)

    if program.lower() == "gaussian":
        mol_char = rf.mol_from_gauss(pop_file, pop=method)
        print(mol_char)
        mol_char.bonding = in_mol.bonding
        mol_char.thresh = in_mol.thresh
        charges = [i.q for i in mol_char]
        print(mol_char.bonding,mol_char.thresh)
        # correct charges if they are not perfectly neutral
        if sum(charges) != 0.0:
            outfile.write("Charge correction: " +
                                str(sum(charges)) + "\n")
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
    def_inputs = {
        "high_level" : "gaussian",
        "low_level"  : "gaussian" }

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

def neutralise_cluster(clust):
    """Remove any net charge from cluster"""

    net_charge = sum([atom.q for atom in clust])
    print("Initial charge: ", net_charge)
    print("Adding correction: ", net_charge/len(clust))
    for atom in clust:
        atom.q -= net_charge/len(clust)
    net_charge = sum([atom.q for atom in clust])
    print("Final charge: ", net_charge)
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


def vis_pce(mol, charges, name):
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
    print("number of charges: ", len(vis))
    vis.remove_duplicates()
    print("number of charges (2): ", len(vis))
    vis.write_xyz(name)
    return


def prep_geoopt(agg, region_2, high_points):
    """
    Hard-coded to gaussian16 and xTB; a more extensible framework can be made
    """
    here = os.getcwd()
    # write fromage directories 
    calc_paths = ["opt", "opt/mh", "opt/ml", "opt/rl"]
    for cpath in calc_paths:
        if not os.path.exists(cpath):
            os.mkdir(cpath)

    high_level_write, low_level_write = getPrograms()

    # write ml
    os.chdir("opt/ml")
    low_level_write("ml.temp", [], region_2, os.path.join(here, "ml.template"))
    os.chdir(here)

    # write rl
    os.chdir("opt/rl")
    low_level_write("rl.temp", region_2, [], os.path.join(here, "rl.template"))
    os.chdir(here)

    # write mh 
    os.chdir("opt/mh")
    high_level_write("mh.temp", [], region_2, os.path.join(here, "mh.template"))
    os.chdir(here)
    return


def run_opt(agg):
    """
    Call Sequence method from fro_run.py and peform geometry optimisation in this code

    
    """



    return

def main(n_agg_states, n_mono_states, nprocs, vis_charges, run_type, parallel_job, custom_region): 
    """
    Run RunSequnce object to get point charges, then create three files: mono1.com, mono2.com and agg

    So far, Ewald is performed for just the agg, then the local charges modified. Indicies is used to specify the fragments for some other scripts I've written
    """
    # location
    print("entering main")
    here = os.getcwd()
    outfile = open(here + "/frooverdia.out", "w")

    # print start time
    start_time = datetime.now()
    outfile.write("STARTING TIME: " + str(start_time) + "\n")

    # read config inputs
    inputs = pcf.parse_inputs("config")

    # read the input cell
    cell = rf.mol_from_file(inputs["cell_file"])
    cell.vectors = inputs["vectors"]
    print(len(cell))
    cell.bonding = inputs["bonding"]
    cell.thresh = inputs["bond_thresh"]
    cell = cell.confined()
    outfile.write("Read " + str(len(cell)) + " atoms in cell_file\n")

    # get indices from config
    indices = inputs["atom_label"]
    if len(indices) < 2:
        raise ValueError("At least two molecules must be included in model for fragment diabatisation")
    n_monomers = len(indices)
    outfile.write(f"Number of fragments: {n_mono_states}")
    print("number of fragments", n_monomers)

    # High level charge assignment to cell
    print("Starting centered_cell")
    populate_cell(cell, inputs["high_pop_program"], inputs["high_pop_file"], inputs["high_pop_method"])
    region_1, cell = cell.centered_mols(indices) ## use fragments specified by argparse
    region_1_pc = region_1.copy()


    region_1.write_xyz("mol.init.xyz")
    if inputs["print_tweak"]:
        ef.write_xyz("tweaked_cell.xyz", cell)


    # check len of region_1 is the same as expected
    if not len(region_1.segregate()) == n_monomers: 
        raise ValueError("Different number of monomers than specified in config")
    
    # generate shell region and PCE
    print("starting RunSeq")
    run_sequence = rs.RunSeq(region_1, cell, inputs)
    region_2, high_points = run_sequence.run()
    region_2.write_xyz("shell.xyz")
    print("Finished RunSeq")

    ## overwrite generated file read-in custom QM region (e.g. defect or optimised)
    print("assinging charges to QM region")
    if custom_region != None:
        region_1 = rf.mol_from_file("model.init.xyz")
        fc.assign_charges(region_1_pc, region_1)


    # neutralise shell and write agg file
    high_points = neutralise_cluster(high_points)

    # divide agg into respective monomers and reconsitute; consistency of atom indexing between (mono1 + mono2) and agg
    # is required for Overdia to run 
    agg_as_mols  = region_1.segregate(diff_mols=False)

    ## iterate to get list of fragments and the aggregate object. NB: the atoms in fragments and monomers *must* be in same order
    monomers, mono_envs, mono_paths, agg = [], [], [], Mol([])
    print(here) 
    
    # iterate through number of fragments
    for i in range(n_monomers):
        
        # get fragments
        monomers.append(agg_as_mols[i])

        # reorder agg to match monomers
        agg.extend(agg_as_mols[i])
            
        # get list of environment cops
        mono_envs.append(neutralise_cluster(high_points.copy()))

        # get file path
        mono_paths.append(f"mono{i+1}.com")

    # Make aggregate gaussian16
    agg_p = "agg.com"
    ef.write_gauss(agg_p, agg, high_points, "pce.temp")
    set_gauss(agg_p, nstates=n_agg_states, type="agg", nprocs=nprocs)

    # visualise
    if vis_charges:
       vis_pce(agg,high_points, "vis/agg.xyz")


    # run geometry optimisation
    if run_type != "sp":
        
        prep_geoopt(agg, high_points, region_2)

        #2. generate fromage calculation input files
        #3. call scipy.optimize.minize

        # rewrite agg_as_mols
        # raise NotImplementedError("Implement me")


    # make the gaussian16 input
    for i,  (mono, mono_env, path) in enumerate(zip(monomers, mono_envs, mono_paths)):
        ## add missing point charges from dimer calculation 
        for j in range(len(agg_as_mols)):
            if j != i:
               print(f"Adding missing charges in mono{i} environment") 
               mono_env.extend(agg_as_mols[j])

        # print("mono1 env: ", len(mono1_env))
        ef.write_gauss(path, mono, mono_env, "pce.temp")
        set_gauss(path, nstates=n_mono_states, type=f"mono{i+1}", nprocs=nprocs)

        # visualise
        if vis_charges:
            vis_pce(mono,mono_env, f"vis/mono{i+1}.xyz")

    outfile.close()

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for generating input files for Overdia diabatisation calculation. Will generate three files (mono1.com, mono2.com, and agg.com) for Gaussian16 calculations. Template file pce.tmp must be specified.")
    
    # Define optional arguments for variable settings
    parser.add_argument('--N', type=int, default=60, help='Make input with N aggregate states in Gaussian16')
    parser.add_argument('--M', type=int, default=10, help='Make input with M monomer states in Gaussian16')
    parser.add_argument('--nprocs', type=int, default=40, help='Request nprocs processors in Gaussian16')
    parser.add_argument('--vis_charges', type=bool, default=True, help='Visualize charges')
    parser.add_argument('--run_type', type=str, default='sp', help='sp: single-point, opt: geometry optimisation')
    parser.add_argument('--parallel_job', type=bool, default=True, help='Run in parallel')
    parser.add_argument('--custom_region', type=str, default=None, help='Custom region file')


    args = parser.parse_args()
    main(args.N, args.M, args.nprocs, args.vis_charges, args.run_type, args.parallel_job, args.custom_region)

                

                

