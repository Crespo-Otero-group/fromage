#!/usr/bin/env python
## Interface between fromag and Newton-X
## 
## Author
## Federico J Hernandez
## October 2023

import time,datetime,os,sys
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
    states = []
    state = None
    natoms = None

    if os.path.isfile("current_state.dat"):
        with open("current_state.dat") as f:
            state = int(f.read().strip())
    

    try:
        with open("../../user_config.nml", "r") as file:
            for line in file:
                line = line.strip()
                if '=' not in line:
                    continue
                key = line.split('=')[0].strip()
                val = line.split('=')[1].strip().rstrip(',')
                if key == 'nat':
                    natoms = int(val)
                elif key == 'nstat':
                    states.append(int(val))
                elif key == 'nstatdyn' and state is None:
                    state = int(val)
    except FileNotFoundError:
        try:
            with open("initqp_input", "r") as data:
                lines = data.readlines()
            for line in lines:
                if "NUMAT" in line:
                    natoms = int(line.split()[2])
                if "NFS" in line:
                    states = [int(line.split()[2])]
                if "NIS" in line:
                    state = int(line.split()[2])
        except FileNotFoundError:
            print("Error in subroutine fromage/utils/newtonx/fro_nx.py.")
            print("Both 'user_config.nml' and 'initqp_input' files not found.")
            print("Check you have properly set the input files for NX")

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

    mol_file = inputs.get("mol_file", "mol.init.xyz")
    if "shell_file" in inputs:
        shell_file = inputs["shell_file"]
    elif "shell_file_flex" in inputs:
        shell_file = inputs["shell_file_flex"]
    else:
        shell_file = "shell.xyz"
    shell_file_fixed = inputs.get("shell_file_fixed", "shell_fixed.xyz")
    high_level = inputs["high_level"]
    low_level = inputs["low_level"]
    if "hl_natoms" in inputs.keys():    
        hl_natoms = int(inputs["hl_natoms"])
    else:
        hl_natoms = None
    if "ll_flex_natoms" in inputs.keys():        
        ll_natoms = int(inputs["ll_flex_natoms"])
    else:
        ll_natoms = None
    nprocs = inputs["nprocs"]

    at_reparam = inputs["at_reparam"]
    if at_reparam:
         at_reparam = []
         at_reparam = [int(x) for x in inputs["at_reparam"]]
         at_reparam = np.array(at_reparam)

    pop_an = inputs["pop_an"]
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
#    nstates = int(np.sum(states))

    for n, s in enumerate(states):
        ms = int(spin[n] * 2 + 1)
        mult.append(ms)

    soc_coupling = []

    vals.append(out_file)
    vals.append(mol_file)
    vals.append(shell_file)
    vals.append(shell_file_fixed)
    vals.append(high_level)
    vals.append(low_level)
    vals.append(hl_natoms)
    vals.append(ll_natoms)
    vals.append(pop_an)
    vals.append(nprocs)
    vals.append(singlestate)
    vals.append(spin)
    vals.append(mult)
    vals.append(soc_coupling)
    vals.append(at_reparam)

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

def write_ener_oos(energies, oos):
    """
    Write energies and oscillator strengths to separate files.
    """
    e_nx  = open("fro_energies.dat", "w")
    oos_nx = open("fro_oos.dat", "w")
    for energy in energies:
        e_nx.write(f"{energy}\n")
    for os in oos:
        oos_nx.write(f"{os}\n")

#    for i in range(1,oos.shape[0]+1):
#        with open("epot.{}".format(i + 1), 'w') as en_fl, \
#             open("oos.{}".format(i + 1), 'w') as oos_fl:
#            en_str = "{:15.10f}\n{:15.10f}".format(
#                energies[0], energies[i])
#            oos_str = "{:15.10f}".format(oos[i-1])
#            en_fl.write(en_str)
#            oos_fl.write(oos_str)
#    subprocess.run("mv mh/epot* .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
#    subprocess.run("mv mh/oos* .", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    return

def write_nx_info(high_level,energies,grads,nacs,socs,oos=[],flex_natoms=None):
    """
    """

    # Write information for initial Conditions
    if len(oos) > 0:
        write_ener_oos(energies,oos)
        return None

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
        nac_ll = 0.
        for i in range(nacs.shape[0]):
            for j in range(nacs.shape[1]):
                nac_str = "{:15.10f} {:15.10f} {:15.10f}".format(
                    nacs[i,j,0], nacs[i,j,1], nacs[i,j,2]) + "\n"
                nx_files[2].write(nac_str)
            if flex_natoms:
                for k in range(flex_natoms):
                    nac_str = "{:15.10f} {:15.10f} {:15.10f}".format(
                        0.0, 0.0, 0.0) + "\n"     
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
    elif method == 'orca':
        copy_outputs_orca()

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

################# COMPLETE HERE ORCA ###################
def copy_outputs_orca():
    """
    """
    subprocess.run("cp mh/mh.out orca.out", shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    return None

########################################################
def _chk_hlevel_in_methods(high_level):
    """
    """
    methods = ['molcas',
               'turbomole',
               'turbomole_tddft',
               'gaussian',
               'orca',
               'fomo-ci',
               'mopac']
    # Check if the high_level method is supported for SH-dynamics with NX
    if high_level in methods:
       pass
    else:
        out_file.write(" The method %s is not implemented in fromage&Newton-X\n" % (high_level))
        out_file.write("The job is dying now :-( ")
        sys.exit("The method %s is not implemented in fromage&Newton-X :-(\n" % (high_level))
    
    return None

def get_mol_shell_atoms(mol_file, shell_file, flex=None, fixed_shell_file=None):
    """
    """
    # read initial coordinates
    mol_atoms = rf.read_xyz(mol_file)[0]

    if flex:
        shell_atoms = rf.read_xyz(shell_file)[0]
    else:
        shell_atoms = []

    fixed_atoms = []
    if flex and fixed_shell_file and os.path.isfile(fixed_shell_file):
        fixed_atoms = rf.read_xyz(fixed_shell_file)[0]

    # make the initial QM coordinates into a flat list
    atoms_array = []
    for atom in mol_atoms:
        atoms_array.append(atom.x)
        atoms_array.append(atom.y)
        atoms_array.append(atom.z)

    in_pos = np.array(atoms_array)

    return in_pos, mol_atoms, shell_atoms, fixed_atoms

def newtonx_initconds(inputs,natoms,states,state):
    """
    """
    (out_file, mol_file, shell_file, shell_file_fixed, high_level,
    low_level, hl_natoms, ll_natoms,  pop_an,
    nprocs, singlestate, spin, mult, soc_coupling,
    at_reparam) = parse_fro_input(inputs,states)

    write_nx_head(out_file)
    _chk_hlevel_in_methods(high_level)

    flex = None
    if hl_natoms and ll_natoms:
        flex = True
        natoms_flex = ll_natoms
        dim_hl = int(3*hl_natoms)

    in_pos, mol_atoms, shell_atoms, fixed_atoms = get_mol_shell_atoms(
        mol_file, shell_file, flex, shell_file_fixed)

    # initialise calculation objects
    rl = calc.setup_calc("rl", low_level)
    ml = calc.setup_calc("ml", low_level)
    mh = calc.setup_calc("mh", high_level)

    pass_nac = []
    in_cond = True
    if flex:
        run_calcs(out_file,in_pos,mh,ml,rl,mol_atoms,at_reparam,pop_an,nprocs,state,states,
                  singlestate,pass_nac,soc_coupling,shell_atoms,hl_natoms,ll_natoms,in_cond,
                  fixed_atoms=fixed_atoms)
    else:
        run_calcs(out_file,in_pos,mh,ml,rl,mol_atoms,at_reparam,pop_an,nprocs,state,states,
                  singlestate,pass_nac,soc_coupling,in_cond=in_cond)
    natoms_mh = hl_natoms if hl_natoms is not None else len(mol_atoms)    
    mh_en_gr = mh.read_out(in_pos,natoms_flex = ll_natoms,natoms = hl_natoms,state = state,
                           states = states, mult = mult, singlestate = singlestate,
                           soc_coupling = soc_coupling, in_cond=in_cond)

    mh_en, mh_gr_tmp, mh_scf, nac, soc = mh_en_gr

    oos = mh.read_osc_str()

    write_nx_info(high_level=high_level,energies=mh_en,grads=[],
                  nacs=[], socs=[], oos=oos)

    return

def newtonx_sequence(inputs,natoms,states,state):
    """
    """
    # Parse fromage.in file
    (out_file, mol_file, shell_file, shell_file_fixed, high_level,
    low_level, hl_natoms, ll_natoms,  pop_an,
    nprocs, singlestate, spin, mult, soc_coupling,
    at_reparam) = parse_fro_input(inputs,states)

    pcgrad_bool = bool_cast(inputs.get("pcgrad", "0"))

    write_nx_head(out_file)
    _chk_hlevel_in_methods(high_level)

    flex = None
    if hl_natoms and ll_natoms:
        flex = True
        natoms_flex = ll_natoms
        dim_hl = int(3*hl_natoms)

    in_pos, mol_atoms, shell_atoms, fixed_atoms = get_mol_shell_atoms(
        mol_file, shell_file, flex, shell_file_fixed)

    print("[DEBUG newtonx_sequence] mol_file={} shell_file={} shell_file_fixed={}".format(
        mol_file, shell_file, shell_file_fixed))
    print("[DEBUG newtonx_sequence] len(mol_atoms)={} len(shell_atoms)={} len(fixed_atoms)={} len(in_pos)={}".format(
        len(mol_atoms), len(shell_atoms), len(fixed_atoms), len(in_pos)))
    print("[DEBUG newtonx_sequence] in_pos[:3]={}".format(in_pos[:3]))

    methods_wnacs = ['molcas', 'dftb'] # Extend this list to other methods that compute NACs
    pass_nac = []

    if high_level in methods_wnacs:
        vdoth = get_vdoth()
        if vdoth == 0:
            pass_nac = get_nacs_coup()

    # initialise calculation objects
    rl = calc.setup_calc("rl", low_level)
    ml = calc.setup_calc("ml", low_level)
    mh = calc.setup_calc("mh", high_level)

    if flex:
        run_calcs(out_file,in_pos,mh,ml,rl,mol_atoms,at_reparam,pop_an,nprocs,state,states,
                  singlestate,pass_nac,soc_coupling,shell_atoms,hl_natoms,ll_natoms,None,
                  fixed_atoms=fixed_atoms)
    else:
        run_calcs(out_file,in_pos,mh,ml,rl,mol_atoms,at_reparam,pop_an,nprocs,state,states,
                  singlestate,pass_nac,soc_coupling)

    # read results. Each x_en_gr is a tuple (energy,gradients,scf_energy)
    if flex:
        rl_en_gr = rl.read_out(in_pos, in_mol=mol_atoms, in_shell=shell_atoms, natoms_flex=ll_natoms)
        if low_level == 'dftb':
            ml_en_gr = ml.read_out(in_pos[:dim_hl], natoms_flex=ll_natoms, pcgrad=pcgrad_bool)
        else:
            ml_en_gr = ml.read_out(in_pos[:dim_hl], natoms_flex=ll_natoms)

        if high_level == 'molcas':
            mh_en_gr = mh.read_out(in_pos, natoms_flex=ll_natoms, natoms=hl_natoms, state=state,
                                   states=states, mult=mult, singlestate=singlestate,
                                   soc_coupling=soc_coupling, newtonx=True)
        else:
            mh_en_gr = mh.read_out(in_pos, natoms_flex=ll_natoms, natoms=hl_natoms, state=state,
                                   states=states, mult=mult, singlestate=singlestate,
                                   soc_coupling=soc_coupling, pcgrad=pcgrad_bool)
        
    else:
        rl_en_gr = rl.read_out(in_pos,in_mol = mol_atoms,in_shell = shell_atoms)
        ml_en_gr = ml.read_out(in_pos)
        if high_level == 'molcas':
            mh_en_gr = mh.read_out(in_pos, natoms_flex = ll_natoms, natoms = natoms, state = state,
                                   states = states, mult = mult, singlestate = singlestate,
                                   soc_coupling = soc_coupling, newtonx = True)
        else:
            mh_en_gr = mh.read_out(in_pos, natoms_flex = ll_natoms, natoms = natoms, state = state,
                                   states = states, mult = mult, singlestate = singlestate,
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
    mh_en, mh_gr_tmp, mh_scf, nac, soc = mh_en_gr
    print("[DEBUG] mh_gr:{}".format(np.array(mh_gr_tmp)))
    ml_en, ml_gr, ml_scf, _, _ = ml_en_gr
    rl_en, rl_gr, rl_scf, _, _ = rl_en_gr
    print("[DEBUG] len(ml_gr)={} len(rl_gr)={}".format(len(ml_gr), len(rl_gr)))
    print("[DEBUG] mh_gr_tmp shape={}".format(np.array(mh_gr_tmp).shape))

    if flex:
        flex_natoms = hl_natoms + ll_natoms
        nstates = int(np.sum(states))
        ml_gr = np.array(ml_gr).reshape((1, flex_natoms, 3))
        mh_gr = np.zeros((nstates, flex_natoms, 3))
        mh_gr[:,:hl_natoms,:] = mh_gr_tmp
        if pcgrad_bool and high_level == "orca":
            pcgrad_file = os.path.join("mh", "mh.pcgrad")
            if os.path.isfile(pcgrad_file):
                pc_grad = rf.read_orca_pcgrad(pcgrad_file)
                pc_grad_flex = pc_grad[:ll_natoms] * bohrconv
                mh_gr[state - 1, hl_natoms:, :] = pc_grad_flex
                print("[pcgrad] mh nuclear grad norm: {:.6e} Ha/Ang  "
                      "shell grad norm: {:.6e} Ha/Ang".format(
                      np.linalg.norm(mh_gr[:, :hl_natoms, :]),
                      np.linalg.norm(pc_grad_flex)))
            else:
                print("[pcgrad] WARNING: mh.pcgrad not found at {}".format(pcgrad_file))
        rl_gr = np.array(rl_gr).reshape((1, flex_natoms, 3))
    else:
        ml_gr = np.array(ml_gr).reshape((1, natoms, 3))
        mh_gr = mh_gr_tmp
        rl_gr = np.array(rl_gr).reshape((1, natoms, 3))
    
    en_combo = rl_en - ml_en + mh_en
    gr_combo = mh_gr.copy()
    gr_combo[state-1,:,:]  = rl_gr - ml_gr + mh_gr[state-1,:,:]
    scf_combo = rl_scf - ml_scf + mh_scf

    if flex:
        write_nx_info(high_level,en_combo,gr_combo,nac,soc,[],ll_natoms)
    else:
        write_nx_info(high_level,en_combo,gr_combo,nac,soc)

    write_ONIOM_info(out_file,mh_en,ml_en,rl_en,en_combo,scf_combo,state)

    return None  

def run_calcs(out_file,in_pos,mh,ml,rl,mol_atoms,at_reparam,pop_an,nprocs,state,states,singlestate,
              pass_nac=[],soc_coupling=[],shell_atoms=None,hl_natoms=None,ll_natoms=None,in_cond=None,
              fixed_atoms=None):
    """
    Run the calculations as subprocesses in parallel
    Note that the schemes are different if frozen or
    flexible ONIOM is required
    """

    calcs = []

    if hl_natoms and ll_natoms:
        all_atoms = mol_atoms + shell_atoms
        dim_hl = 3 * hl_natoms

        flex_pos_array = []
        for atom in shell_atoms:
            flex_pos_array.extend([atom.x, atom.y, atom.z])
        all_pos = np.concatenate((in_pos, np.array(flex_pos_array)))

        rl_proc = rl.run(atoms=ao.array2atom(all_atoms, all_pos), nprocs=nprocs)
        rl_proc.wait()
        rl_charges_array = rl.read_charges(pop=pop_an)

        flex_pc = ao.array2atom(shell_atoms, np.array(flex_pos_array),
                                rl_charges_array[hl_natoms:hl_natoms + ll_natoms])
        if fixed_atoms:
            fixed_pos_array = []
            for atom in fixed_atoms:
                fixed_pos_array.extend([atom.x, atom.y, atom.z])
            fixed_pc = ao.array2atom(fixed_atoms, np.array(fixed_pos_array),
                                     rl_charges_array[hl_natoms + ll_natoms:
                                                      hl_natoms + ll_natoms + len(fixed_atoms)])
            all_pc = flex_pc + fixed_pc
        else:
            all_pc = flex_pc

        mh_proc = mh.run(ao.array2atom(mol_atoms, in_pos[:dim_hl]),
                         all_pc,
                         nprocs,
                         state=state,
                         states=states,
                         singlestate=singlestate,
                         nac_coupling=pass_nac,
                         soc_coupling=soc_coupling)
        calcs.append(mh_proc)
        mh_proc.wait()
        if in_cond is None:
            ml_proc = ml.run(ao.array2atom(mol_atoms, in_pos[:dim_hl]),
                             all_pc,
                             nprocs)
            calcs.append(ml_proc)
            ml_proc.wait()
    else:
        mh_proc = mh.run(atoms=ao.array2atom(mol_atoms, in_pos),
                         nprocs=nprocs,
                         state=state,
                         states=states,
                         singlestate=singlestate,
                         nac_coupling=pass_nac,
                         soc_coupling=soc_coupling)
        calcs.append(mh_proc)
        if in_cond is None:
            rl_proc = rl.run(atoms = ao.array2atom(mol_atoms, in_pos), nprocs = nprocs)
            calcs.append(rl_proc)
            ml_proc = ml.run(atoms = ao.array2atom(mol_atoms, in_pos),nprocs = nprocs)
            calcs.append(ml_proc)

        # Wait until all parallel calculations are finished
        for proc in calcs:
            proc.communicate()

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

