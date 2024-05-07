#!/usr/bin/env python
"""
Implementation of the EE(L1) and EE(L2) point-charge embedding scheme of Caricato and
co-workers (https://doi.org/10.1063/1.4972000). 

This scheme uses constrained ESP fitting to improve the envionrment point charges ({qR} in their paper's 
vernacular).
Following an initial rl calculation, the charges of the central region {qM} are kept fixed, and the LACs {qB} 
which connect atoms {qR} and {qM} are extracted. Next, mh and ml calculations are performed, from which the link
atoms provide new charges at the high {qLinkH} and low {qLinkL} levels. These are used to update {qB} by taking the 
average as in EE(L1), or seperately, as in EE(L2). A constrained RESP fit then optimises the {qR charges}. This is 
repeated until {qR} charges are self-consistent. 

I'm keeping this as a seperate script to my other link atom point charge schemes, as it requires
running Gaussian/psi4. I think I'm going to implement this using CP2K

Initial development 11/08/23
Author: Michael Ingham

This is going to require a lot of use of subprocess

psuedo-code:

1.) Run rl                                      # for now this will be gaussian
2.) get rl.fchk                                 # inp for multiwfn
3.) resp fit rl using multiwfn                  # get charges
4.) partition charges into {qM}, {qB}, {qR}     # charge by region
5.) run ml, mh embedded in {qR}                 
6.) extract {qLinkL} from ml
7.) extract {qLinkH} from mh
8.) {qB} = 1/2({qLinkL}+{qLinkH})               # EE(L1), treat seperately for EE(L2)
9.) generate plain text file        # this is slightly complicated, it will need to read all the charges 
                                    # from {qM} and {qB} as constraints to their correspodning position in 
                                    # the calculation, the cluster must also be constrained to beneutral/integer 
                                    # charge
10.) contrained RESP on rl using multwifn fixing 
11.) obtain new {qR}, {q'R}
12.) check if {qR} = {q'R}
13.) yes: stop; no: go back to 5.)

"""
import fromage.io.edit_file as ef
import fromage.io.read_file as rf
import fromage.scripts.fro_assign_charges as fs 
from fromage.io import parse_config_file as pcf
from fromage.utils.linkatom import LinkAtom
from fromage.utils.mol import Mol
from fromage.utils import calc
import numpy as np
import subprocess
import argparse
import os 
from datetime import datetime
from statistics import mean
import fromage.utils.calc as calc
import fromage.utils.array_operations as ao


def make_mwfn_constraint(qM_charges, qB_charges):
    """makes constraints file for RESP fit"""
    constraint_file_path = "constraint_file"

    #remove previous
    if os.path.exists(constraint_file_path):
        os.remove(constraint_file_path)

    with open(constraint_file_path, "w") as constraint_file:
        for i, char in enumerate(qM_charges): # this stays fixed
            constraint_file.write(f"{i+1}, {char.q}\n")  # atom number, constraint value
        for i, char in enumerate(qB_charges): # this changes every iteration
            constraint_file.write(f"{i+len(qM_charges)+1}, {char.q}\n")
    return


### USEFUL BITS

    # real_file = "rl.chg"
    # real_file_in = "rl.molden"

    # nprocs = 32
    # MK_grid = 1.0 # 1 is likely too low; this parameter massively influences the efficiency of the
    #               # of the multiwfn RESP calculation, default is 6.0. See QMOF-d29cec2 benchmark  
    #               # RESP using ChelpG-fit  
    # auto_constrained = f"7\n18\n3\n1\n1\n{MK_grid}\n0\n6\n1\nconstraint_file\n2\n1.22\ny\n0\n0\nq"
    # auto_resp_file = f"7\n18\n3\n1\n1\n{MK_grid}\n0\n1\n1.22\n\ny\n0\n0\nq"

   #get input file for multiwfn
    # make_mwfn_config(name="auto_resp", mwfn_inp=auto_resp_file)


def main(mol, real, program="xtb", scheme="L2"):
    """redistribute charges according to the EE(L1) or EE(L2) charge redistribution schemes"""

    # initialise calculation
    here = os.getcwd()


    with open(here+"output", "w") as output_file:
        # print start time
        start_time = datetime.now()
        output_file.write("STARTING TIME: " + str(start_time))

    if os.path.exists("RESP.dat"):
        os.remove("RESP.dat")

    # useful bits
    if scheme=="L2":
        low_level= None
        high_level = args.program
    
    # tolerance
    q_thresh = 1e-10

    #get initial files
    mol = rf.mol_from_file(mol)
    real = rf.mol_from_file(real)

    mol.set_bonding(args.bonding, args.thresh)
    real.set_bonding(args.bonding, args.thresh)
    real_array = []
    for atom in real:
        real_array.append(atom.x)
        real_array.append(atom.y)
        real_array.append(atom.z)
    real_array = np.array(real_array)

    # get shell region
    shell = real.remove_model(mol)
    shell.set_bonding(args.bonding, args.thresh)
    shell.write_xyz("shell.xyz")

    # make aug model
    atoms_lac, atoms_lah = shell.detect_bondcuts(mol)

    mol = mol.rearrange_mol(atoms_lac)
    real = Mol(mol+ atoms_lah +shell) # redefine real region in useful order 
    real.remove_duplicates()
    print(len(real))
    aug_mol, linkatoms = mol.add_linkatoms(atoms_lac, atoms_lah)
    aug_mol.write_xyz("aug_mol.init.xyz")


    # make directory for rl calculation
    rl_path = "rl_ee"
    if not os.path.exists(rl_path):
        os.makedirs(rl_path)

    #run xtb calculation to get real charges
    rl_ee = calc.xtb_calc("rl_ee")
    rl_ee.run_quick(real,"rl_ee")

    # Run Multiwfn_no_GUI command with auto_resp as input
    print("Starting rl-RESP fit")
    rl_ee_mwfn = calc.Multiwfn_calc("molden.input", "rl_ee")
    rl_ee_mwfn.set_input(type="RESP")
    rl_ee_mwfn.run(args.nprocs)
    rl_resp = rl_ee_mwfn.read_charges()

    print("initial resp charges:", rl_resp)
    #assign RESP charges to real region (assumes input in same order)
    for i, atom in enumerate(real): 
        atom.q = rl_resp[i]




    iteration = 0
    scf_active = True
    while scf_active:
        
        iteration +=1
        print(f"\n------------------\nIteration {iteration}\n------------------\n")

            # #get charges
        qB_len = len(atoms_lac)
        qM_len = len(mol)
        qR_len = len(real) - qM_len - qB_len

        qM_charges = real[:qM_len]                    # charges on model region atoms
        qB_charges = real[qM_len:qM_len+qB_len]   # charges on LAH 
        qR_charges = real[qM_len+qB_len:]             # charges on remaining atoms in system

        print("{qM}", len(qM_charges), "{qB}", len(qB_charges), "{qR}", len(qR_charges), "total", len(qM_charges)+len(qB_charges)+len(qR_charges) ) # keep fixed
        #print("qM initial: ",qM_charges)
        print("qB initial: ",qB_charges)
        # # create directories for model calculat   ions
        # mh_path = os.path.join(here, 'mh') 
        # if not os.path.exists(mh_path):
        #     os.makedirs(mh_path)

        # model low calculation
        ml_path = os.path.join(here, 'ml_ee') 
        if not os.path.exists(ml_path):
            os.makedirs(ml_path)
        ml_ee = calc.xtb_calc("ml_ee")
        ml_ee.run_quick(aug_mol,"ml_ee", charges=qR_charges)

        # model high calculation
        mh_path = os.path.join(here, 'mh_ee') 
        if not os.path.exists(mh_path):
            os.makedirs(mh_path)
        mh_ee = calc.Gauss_calc("mh_ee")

        mh_ee.write_g_temp("mh.temp", [], qR_charges, os.path.join(here, "mh.template"))
        mh = calc.setup_calc("mh", high_level)
      


        # get new charges
        ml_ee_mwfn = calc.Multiwfn_calc("molden.input", "ml_ee")
        ml_ee_mwfn.set_input(type="RESP")
        ml_ee_mwfn.run(args.nprocs)
        ml_charges = ml_ee_mwfn.read_charges()

        #get store for now
        #mh_charges = ml_charges.copy()

        # print("Starting model-high calculation")
        # os.chdir(mh_path)
        # mh_temp_path = os.path.join(mh_path, "mh.temp")
        # mh_com_path = os.path.join(mh_path, "mh.com")

        # if os.path.exists(mh_temp_path):
        #     os.remove(mh_temp_path)

        # if os.path.exists(mh_com_path):
        #     os.remove(mh_com_path)
         
        # high_level_write("mh.temp", [], qR_charges, os.path.join(here, "mh.template"))
        # mh = calc.setup_calc("mh", high_level)
        # os.chdir(here)
    
        #run model calculations
        #print("Running initial gaussian calculation")
        # if low_level == "gaussian":
        #     ml_proc = ml.run(mol_aug_model, nprocs)
        #     ml_proc.wait()
        #     mh_proc = mh.run(mol_aug_model, nprocs)
        #     mh_proc.wait()
        # elif low_level =="xtb":
        #     print("Running xTB calculation")
        #     xtb_command = "xtb aug_model.xyz --scc --molden > xtb.out"
        #     os.chdir("./ml")
        #     xtb_proc = subprocess.run(xtb_command, shell=True, check=True)
        #     #xtb_proc.wait()
        #     os.chdir(here)
        #     print("Done xTB caclulation")
        #     print("Running Gaussian calculation")
        #     mh_proc = mh.run(mol_aug_model, nprocs)
        #     mh_proc.wait()
        #     print("Done Gaussian")


        # mh RESP fit
        # print("Obtaining formatted checkpoint files for mh calculation:")
        # if high_level == "gaussian":
        #     subprocess.call(['formchk', 'mh/mh.chk'])
        #     mwfn_run("auto_resp", "mh/mh.fchk", nprocs)
        # else:
        #     raise NotImplementedError(f"model-high RESP fit not implemented in {high_level}") 

        # # ml RESP fit
        # print("Obtaining formatted checkpoint files for ml calculation:")
        # if low_level == "gaussian":
        #     subprocess.call(['formchk', 'ml/ml.chk'])
        #     mwfn_run("auto_resp", "ml/ml.fchk", nprocs)
        # elif low_level == "xtb":
        #     subprocess.run(["mv", "ml/molden.input", "ml/ml.molden"])
        #     mwfn_run("auto_resp", "ml/ml.molden", nprocs)
        
        # get link atom charges
        qLinkH_charges = mh_charges[-len(qB_charges):]
        qLinkL_charges = ml_charges[-len(qB_charges):]
        print("{qLinkH}", qLinkH_charges)
        print("{qLinkL}", qLinkL_charges)

        # update qB as per EE(L1) scheme
        averaged_charges = [
            0.5*(qLinkH+qLinkL) for qLinkH,qLinkL in zip(qLinkH_charges,qLinkL_charges) 
        ]

        print("Averaged charges (EE(L1) scheme):", averaged_charges)
        for i, char in enumerate(qB_charges):
            print(i, char)
            char.q = averaged_charges[i]
        
        #constrained RESP fit
        rl_ee_mwfn = calc.Multiwfn_calc("molden.input", "rl_ee", qM=qM_charges, qB=qB_charges)
        rl_ee_mwfn.set_input(type="RESP")
        rl_ee_mwfn.run(args.nprocs)
        print(rl_resp)
        rl_resp = rl_ee_mwfn.read_charges()
        print(rl_resp)

      
        # update qR charges 
        real_old = real.copy()
        for atom, char in zip(real, rl_resp):
            atom.q = char


        dQ = [abs(old_char.q-new_char.q) for old_char, new_char in zip(real_old, real)]

        print(rl_resp[len(mol):])
        print("Old charges:", [round(char.q,10) for char in real_old][len(mol):])
        print("New charges:",  [round(char.q,10) for char in real][len(mol):])
        print("dQ: ", [round(q,4) for q in dQ][len(mol):])

        dQ = dQ[len(mol):]
        ave_deltaq = mean(dQ)
        max_deltaq = max(dQ)
        
        print("Largest dQ: ", max_deltaq)
        print("Average dQ:", ave_deltaq)
        
        with open("EE_L1.out", "a") as output_file:
                output_file.write("\n-----------------------------") 
                output_file.write("\nIteration: " + str(iteration))
                output_file.write("\nAverage dQ: "  + str(round(ave_deltaq,5)))
                output_file.write("\nLargest dQ: "  + str(round(max_deltaq,5)))
        with open("RESP.dat", "a") as data_file:
            dQ = [round(q,4) for q in dQ]
            data_file.write(f"Iteration {iteration}" + ','.join(map(str, dQ)) + '\n')
        #check converged
        if  ave_deltaq <= q_thresh:
            scf_active = False
            with open("EE_L1.out",  "a") as output_file:
                output_file.write("\n-----------------------------") 
                output_file.write("\nSCF CONVERGED")
                end_time = datetime.now()
                output_file.write("\nELAPSED TIME: " + str(end_time - start_time))
                output_file.write("\nENDING TIME: " + str(end_time))

            print("SCF converged")
        else:
            print("Not converged, next iteration")

        
if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Your program description here.")

    parser.add_argument("--mol", default="mol.init.xyz", help="Path to the model region XYZ file.")
    parser.add_argument("--real", default="real.xyz", help="Path to the real region XYZ file.")
    parser.add_argument("--program", default="gaussian", help="High-level program.")
    # parser.add_argument("--real_chg", default="rl.chg", help="Path to the real charge file.")
    # parser.add_argument("--real_molden", default="rl.molden", help="Path to the real molden file.")
    parser.add_argument("--bonding", default="dis", help="Bonding")
    parser.add_argument("--thresh", type=float, default=1.8, help="Bonding threshold")
    parser.add_argument("--nprocs", type=int, default=12, help="Number of processors to use")
    parser.add_argument("--mk_grid", type=float, default=1.0, help="MK grid parameter.")

    args = parser.parse_args()
    main(args.mol, args.real)
