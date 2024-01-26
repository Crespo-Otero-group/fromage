#!/usr/bin/env python
## Interface between fromag and Newton-X
## 
## Author
## Federico J Hernandez
## October 2023

import time,datetime,os
import numpy as np
import subprocess
import configparser
from fromage.dynamics.periodic_table import Element
from fromage.utils import calc
from fromage.utils import array_operations as ao
from fromage.io import read_file as rf
from fromage.utils.atom import Atom
from fromage.io.parse_config_file import bool_cast

bohrconv = 1.88973  # Something in Angstrom * bohrconv = Something in Bohr

def read_nx_control():
    """
    This subrutine identifies if the Newton-X calculation is to compute a spectra
    or dynamics
    """
    try:
        with open("control.d", "r") as file:
            line = file.readline()
            info = line.split(',')
            natoms = int(info[0].strip())
            states_tmp = info[2].strip()
            states = [int(x) for x in states_tmp]
            nstates = int(np.sum(states))
            state = int(info[3].strip())        
    except FileNotFoundError:
        try:
            config = configparser.ConfigParser()
            config.read("initqp")

            if "dat" in config:
                section = config["dat"]
                if "NUMAT" in section:
                    natoms = int(section["NUMAT"])
                if "NFS" in section:
                    states_tmp = section["NFS"]
                    states = [int(x) for x in states_tmp]
                    nstates = int(np.sum(states))
                if "NIS" in section:
                    state = int(section["NIS"])
             
        except FileNotFoundError:
            print("Error in subroutine fromage/utils/newtonx/fro_nx.py.")
            print("Both 'control.d' and 'initqp' files not found.")
            print("Check you have properly set the input files for NX")
            data = None

    return natoms, states, state

def get_vdoth():
    """
     This subroutined read the sh.inp file to identiy wether the NACs have
     have to computed by the electronic structure spftware
    """
    try:
        with open("sh.inp") as data:
            lines = data.readlines()

    except FileNotFoundError:
        print("Error in subroutine fromage/utils/newtonx/fro_nx.py.")
        print("sh.inp file not found.")

    for line in lines:
        if "vdoth" in line:
            vdoth = int(line.split()[-1])

    return vdoth
    
def get_nacs_coup():
    """
    Read NACs pair from the NX file transmomin to pass them to the software
    that will compute them
    """
    try:
        with open("transmomin") as data:
            lines = data.readlines()

    except:
            print("transmomin file not found - No NACs will be computed.")

    nac_coupling = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            st2 = int(parts[1]) - 1
            st1 = int(parts[3]) - 1
            nac_coupling.append([st2, st1])

    return nac_coupling

def parse_fro_input(inputs,states):

    vals = []
    out_file = inputs["out_file"]
    out_file = open(out_file,"a+")

    mol_file = inputs["mol_file"]
    shell_file = inputs["shell_file"]
    high_level = inputs["high_level"]
    low_level = inputs["low_level"]
    nprocs = inputs["nprocs"]
    if "singlestate" in inputs.keys():
        singlestate = int(inputs["singlestate"])
    else:
        singlestate = 0

    spin = [0]
    statemult = []

    if "spin" in inputs.keys():
        spin = [int(x) for x in inputs["spin"]]

    mult = []
    nstates = int(np.sum(states))

    for n, s in enumerate(states):
        ms = int(spin[n] * 2 + 1)
        mult.append(ms)

    vals.append(out_file)
    vals.append(mol_file)
    vals.append(shell_file)
    vals.append(high_level)
    vals.append(low_level)
    vals.append(nprocs)
    vals.append(singlestate)
    vals.append(spin)
    vals.append(mult)
    vals.append(soc_coupling)

    return (vals)

def write_nx_head(out_file):
    out_file.write("     ###################################\n")
    out_file.write("     #  The Newton-X option is active  #\n")
    out_file.write("     #  fromage is therfore used as a  #\n")
    out_file.write("     # third-party program of Newton-X #\n")
    out_file.write("     ###################################\n")
    return None

def open_files_for_nx(write_nac=False,write_soc=False):
    """
    """
    nx_files = []
    e_nx  = open("fro_energies.dat", "w")
    nx_files.append(e_nx)
    gr_nx = open("fro_gradients.dat", "w")
    nx_files.append(gr_nx)
    if write_nac:
        nac_nx = open("fro_nacs.dat", "w")
        nx_files.append(nac_nx)
    if write_soc:
        soc_nx = open("fro_socs.dat", "w")
        nx_files.append(soc_nx)

    return nx_files

def write_nx_info(high_level,energies,grads,nacs,socs):
    """
    """
    write_nac = False
    write_soc = False
    # Open output files
    if len(nacs.shape) > 2:
        write_nac = True
    if len(socs.shape) > 2:
        write_soc = True
    nx_files = open_files_for_nx(write_nac,write_soc)
    for energy in energies:
        nx_files[0].write(f"{energy}\n")
    
    # Write gradients in Hartree/Bohr 
    grads /= bohrconv
    for i in range(grads.shape[0]):
        for j in range(grads.shape[1]):
            grad_str = "{:15.10f} {:15.10f} {:15.10f}".format(
                grads[i,j,0], grads[i,j,1], grads[i,j,2]) + "\n"
            nx_files[1].write(grad_str)

    if write_nac:
        for i in range(nacs.shape[0]):
            for j in range(nacs.shape[1]):
                nac_str = "{:15.10f} {:15.10f} {:15.10f}".format(
                    nacs[i,j,0], nacs[i,j,1], nacs[i,j,2]) + "\n"
                nx_files[2].write(nac_str)

    # Organise outputs for NX
    copy_output_files(high_level)

    return None
 
def copy_output_files(method):
    """
    This function organise and copy the output files required for 
    NX to do postrun analysis (like cioverlap) outside the mh/ dir
    """
    if method == 'gaussian':
        copy_outputs_gau()
    elif method == 'turbomole' or method == 'turbomole_tddft':
        subprocess.run("cp mh/* .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    elif method == 'molcas':
        copy_outputs_molcas()

    return None

def copy_outputs_gau():
    """
    """    
    subprocess.run("cp mh/*.chk gaussian.chk", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    subprocess.run("mv mh/*.rwf gaussian.rwf", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    subprocess.run("cp mh/mh.com gaussian.com", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    subprocess.run("cp mh/mh.log gaussian.log", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)

    return None

def copy_outputs_molcas():
    """
    """
    subprocess.run("cp mh/molcas/*.log .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    subprocess.run("cp mh/molcas/*.rasscf.molden .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    subprocess.run("cp mh/molcas/*.RasOrb .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        subprocess.run("cp mh/molcas/*.JobIph .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    except:
        print("fromage did not find the .JobIph file")

    return None

def newtonx_sequence(atoms_array,inputs,natoms,states,state):
    """
    """
    # Parse fromage.in file
    (
    out_file,
    mol_file,
    shell_file,
    high_level,
    low_level,
    nprocs,
    singlestate,
    spin,
    mult,
    soc_coupling
    ) = parse_fro_input(inputs,states)

    write_nx_head(out_file)

    methods = ['molcas','turbomole','turbomole_tddft','gaussian']
    methods_wnacs = ['molcas'] # Extend this list to other methods that compute NACs  
  
    # Check the method selected for the high_level is supported for SH-dynamics
    if high_level in methods:
       pass 
    else:
        out_file.write(" The method %s is not implemented for fromage with Newton-X\n" % (high_level))
        out_file.write("The job is dying now :-( ")
        import sys
        sys.exit()

    # read initial coordniates
    mol_atoms = rf.read_xyz(mol_file)[0]

    # read shell atoms
    shell_atoms = rf.read_xyz(shell_file)[0]
    # make the initial coordinates into a flat list
    atoms_array = []
    for atom in mol_atoms:
        atoms_array.append(atom.x)
        atoms_array.append(atom.y)
        atoms_array.append(atom.z)
  
    in_pos = np.array(atoms_array)

    natoms_flex = None
    pass_nac = []

    if high_level in methods_wnacs:
        vdoth = get_vdoth()
        if vdoth == 0:
            pass_nac = get_nacs_coup()


    # initialise calculation objects
    rl = calc.setup_calc("rl", low_level)
    ml = calc.setup_calc("ml", low_level)
    mh = calc.setup_calc("mh", high_level)
 
    # Run the calculations as subprocesses in parallel
    calcs = []    

    if high_level == "fomo-ci" or high_level == "mopac" and at_reparam is not None:
        mh_proc = mh.run(atoms = ao.array2atom(mol_atoms, in_pos),
                         nprocs = nprocs, at_reparam = at_reparam)
    else:
        mh_proc = mh.run(atoms=ao.array2atom(mol_atoms, in_pos), 
                         nprocs = nprocs, 
                         state = state,
                         states = states,
                         singlestate = singlestate,
                         nac_coupling = pass_nac, 
                         soc_coupling = soc_coupling)
    calcs.append(mh_proc)
    if low_level == "fomo-ci" or low_level == "mopac" and at_reparam is not None:
        rl_proc = rl.run(atoms = ao.array2atom(mol_atoms, in_pos), nprocs = nprocs ,at_reparam = at_reparam)
        calcs.append(rl_proc)
        ml_proc = ml.run(atoms = ao.array2atom(mol_atoms, in_pos),
                         nprocs = nprocs , at_reparam = at_reparam)
        calcs.append(ml_proc)
    else:
        rl_proc = rl.run(atoms = ao.array2atom(mol_atoms, in_pos), nprocs = nprocs)
        calcs.append(rl_proc)
        ml_proc = ml.run(atoms = ao.array2atom(mol_atoms, in_pos),nprocs = nprocs)
        calcs.append(ml_proc)

    # Wait until all parallel calculations are finished
    for proc in calcs:
        proc.communicate()

    # read results. Each x_en_gr is a tuple (energy,gradients,scf_energy)
    rl_en_gr = rl.read_out(in_pos,in_mol = mol_atoms,in_shell = shell_atoms)
    ml_en_gr = ml.read_out(in_pos)
#
    if high_level in methods:
        mh_en_gr = mh.read_out(in_pos,
                               natoms_flex = natoms_flex, # FJH
                               natoms = natoms,
                               state = state,
                               states = states,
                               mult = mult,
                               singlestate = singlestate,
                               soc_coupling = soc_coupling)

    """ data format
    mh_en_gr
        0 energy    (nstate, )
        1 grad      (nstate, natoms, 3)
        2 gr_energy float
        3 nac       (nnac, natoms, 3)
        4 soc       (nsoc,)

    ml_en_gr
        0 energy    float
        1 grad      (natoms * 3,)
        2 gr_energy float

    rl_en_gr      
        0 energy    float
        1 grad      (natoms * 3,)
        2 gr_energy float
    """
    mh_en, mh_gr, mh_scf, nac, soc = mh_en_gr
    ml_en, ml_gr, ml_scf = ml_en_gr
    rl_en, rl_gr, rl_scf = rl_en_gr

    ml_gr = np.array(ml_gr).reshape((1, natoms, 3))
    rl_gr = np.array(rl_gr).reshape((1, natoms, 3))

    en_combo = rl_en - ml_en + mh_en
    gr_combo = rl_gr - ml_gr + mh_gr
    scf_combo = rl_scf - ml_scf + mh_scf

    write_nx_info(high_level,en_combo,gr_combo,nac,soc)

    write_ONIOM_info(out_file,mh_en,ml_en,rl_en,en_combo,scf_combo,state)

    return None  

def write_ONIOM_info(out_file,mh_en,ml_en,rl_en,en_combo,scf_combo,state):
    """
    """
    evconv = 27.2114 # eV / Hartree
    # print some updates in the output
    out_file.write("------------------------------\n")
    out_file.write("Real low energy: {:>30.8f} eV{:>14.8f} au\n".format(
        rl_en*evconv, rl_en))
    out_file.write("Model low energy: {:>29.8f} eV{:>14.8f} au\n".format(
        ml_en*evconv, ml_en))
    out_file.write("Model high energy: {:>28.8f} eV{:>14.8f} au\n".format(
        float(mh_en[state-1])*evconv, float(mh_en[state-1])))
    out_file.write(
        "ONIOM Total energy: {:>27.8f} eV{:>14.8f} au\n".format(
        float(en_combo[state-1])*evconv, float(en_combo[state-1])))
    out_file.write(
        "ONIOM SCF energy: {:>29.8f} eV{:>14.8f} au\n".format(
        float(scf_combo) * evconv, float(scf_combo)))
    out_file.write("Gap: {:>42.8f} eV{:>14.8f} au\n".format(
            (float(en_combo[state-1]) - float(scf_combo))*evconv,
            float(en_combo[state-1]) - float(scf_combo)))

    out_file.flush()    

    return None
