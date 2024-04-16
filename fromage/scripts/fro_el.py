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
import decimal
import fromage.io.edit_file as ef
import fromage.io.read_file as rf
import fromage.scripts.fro_assign_charges as fs 
from fromage.io import parse_config_file as pcf
from fromage.utils.linkatom import LinkAtom
from fromage.utils.mol import Mol
from fromage.utils import calc
import numpy as np
import subprocess
import time
import os 
from datetime import datetime
from statistics import mean


def make_mwfn_config(name, mwfn_inp):
    """make auto_resp file for multiwfn calculation"""
    with open(name, "w") as file:
        file.write(mwfn_inp)
    return

def mwfn_run(inp_file, output_file, nprocs, shell_command="Multiwfn_noGUI"):
    """Runs multwfn for given input file"""
    os.environ["np"] = nprocs
    mfwn_command = [shell_command, output_file, f"-nt $np"]
    with open(inp_file, "r") as input_file:
        subprocess.run(mfwn_command, stdin=input_file)
    return

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

def detect_dangling(self):
    """
    Identifies link atom connects (LACs) and link atom hosts (LAHs) by identifying dangling bonds 
    in the molecule. Returns these as seperate molecules so they can be written as input files for fro_poly_run.py

    **CAUTION**: cannot be used for systems with C=O bonds, or unusual bridged systems like B2H6
    
    Parameters
    ----------
    self: Mol
        The molecule containing link atom hosts from which the LACs 
    
    """
    from fromage.utils.mol import Mol
    term_types = ["H", "Cl", "F", "Br", "O"] # etc..
    lah_out = Mol([])
    lac_out = Mol([])
    for atom_a in self:
        tmp_count = 0
        for atom_b in self:
            if self.bonded(atom_a, atom_b) and atom_a != atom_b:
                tmp_count +=1
        if tmp_count == 1 and atom_a.elem not in term_types: # obvs needs to be more general than this
            lah_out.append(atom_a)
            for atom_b in self:
                if self.bonded(atom_a, atom_b) and atom_a != atom_b:
                    lac_out.append(atom_b)

    return lah_out, lac_out

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

def rearrange_model_region(mol_model, atoms_lac,atoms_lah):
    """rearrange mol object"""
    mol_rearranged = Mol([])
    for atom in mol_model:
        if atom not in atoms_lac and atom not in atoms_lah:
            mol_rearranged.append(atom)
    for atom in atoms_lac:
        mol_rearranged.append(atom)
    for atom in atoms_lah:
        mol_rearranged.append(atom)   
    mol_model = mol_rearranged.copy()
    print("Model region:" , len(mol_model), type(mol_model))
    return mol_model

def make_aug_model(mol_model, atoms_lac,atoms_lah):
    """make augmented model region from model and link atoms connects/hosts"""
    mol_aug_model = mol_model.copy()
    print("Making augmented model region", type(mol_aug_model))
    for lac, lah in zip(atoms_lac, atoms_lah):
        link_atom = LinkAtom(lac_partner=lac, lah_partner=lah)
        link_atom.initial_pos()
        print("LAC: ", lac)
        print("LAH: ", lah)
        print("New link atom: ", link_atom)
        if link_atom not in mol_aug_model:
            mol_aug_model.remove(lah)
            mol_aug_model.append(link_atom)
    mol_aug_model.write_xyz("aug_model.xyz")
    return mol_aug_model
def main():
    #inputs = pcf.parse_inputs("config")

    #get cwd
    here = os.getcwd()
    with open(here+"/scc_ee.out", "w") as output_file:
        # print start time
        start_time = datetime.now()
        output_file.write("STARTING TIME: " + str(start_time))

    # generic write file object
    high_level_write, low_level_write = getPrograms()

    #get initial files
    mol_real = rf.mol_from_file("real.xyz")
    mol_model = rf.mol_from_file("model.xyz")
    high_level = "gaussian"
    low_level = "xtb"

    real_file = "rl.chg"
    real_file_in = "rl.molden"

    nprocs = 32
    MK_grid = 1.0 # 1 is likely too low; this parameter massively influences the efficiency of the
                  # of the multiwfn RESP calculation, default is 6.0. See QMOF-d29cec2 benchmark  
                  # RESP using ChelpG-fit  
    auto_constrained = f"7\n18\n3\n1\n1\n{MK_grid}\n0\n6\n1\nconstraint_file\n2\n1.22\ny\n0\n0\nq"
    auto_resp_file = f"7\n18\n3\n1\n1\n{MK_grid}\n0\n1\n1.22\n\ny\n0\n0\nq"

    #get link atoms 
    atoms_lah, atoms_lac = detect_dangling(mol_model)
    atoms_lah.write_xyz("lah.xyz")
    atoms_lac.write_xyz("lac.xyz")
    
    #Rearrange model file so that LACs/LAHs at end
    mol_model = rearrange_model_region(mol_model,atoms_lac,atoms_lah)

    #make augmented model region (with link atoms instead of LAH)
    mol_aug_model = make_aug_model(mol_model,atoms_lac,atoms_lah)

    #get input file for multiwfn
    make_mwfn_config(name="auto_resp", mwfn_inp=auto_resp_file)

    if not os.path.isfile(real_file):
        print("Performing initial RL RESP-fit")
        mwfn_run("auto_resp", real_file_in, nprocs)

    # Run Multiwfn_no_GUI command with auto_resp as input
    print("RL charges obtained for rl")
    rl_resp_charges = [float(line.split()[4]) for line in open(real_file) if line.strip()]   
    print("RL resp charges: ", rl_resp_charges)
    print("Total RESP charge: ", sum(rl_resp_charges))

    #assign RESP charges to real region (assumes input in same order)
    for i, atom in enumerate(mol_real): 
        atom.q = rl_resp_charges[i]

    # #get charges
    qB_len = len(atoms_lac)
    qM_len = len(mol_model) - qB_len
    qR_len = len(mol_real) - qM_len - qB_len

    qM_charges = mol_real[:qM_len]
    qB_charges = mol_real[qM_len:qM_len+qB_len]
    qR_charges = mol_real[qM_len+qB_len:]

    print("{qM}", len(qM_charges), "{qB}", len(qB_charges), "{qR}", len(qR_charges) ) # keep fixed

    # create directories for model calculations
    mh_path = os.path.join(here, 'mh') 
    if not os.path.exists(mh_path):
        os.makedirs(mh_path)

    ml_path = os.path.join(here, 'ml') 
    if not os.path.exists(ml_path):
        os.makedirs(ml_path)

    q_thresh = 1e-4
    #biggest_value = 1000
    iteration = 0
    scf_active = True
    while scf_active:
        iteration +=1
        print(f"\n------------------\nIteration {iteration}\n------------------\n")

        # prepare calculation files
        print("Preparing model-low calculation")
        os.chdir(ml_path)

        ml_temp_path = os.path.join(ml_path, "ml.temp")
        ml_com_path = os.path.join(ml_path, "ml.com")
        


        if os.path.exists(ml_temp_path):
            os.remove(ml_temp_path)

        if os.path.exists(ml_com_path):
            os.remove(ml_com_path)
         
        low_level_write("ml.temp", [], qR_charges, os.path.join(here, "ml.template"))
        ml = calc.setup_calc("ml", low_level)


        print("Starting model-high calculation")
        os.chdir(mh_path)
        mh_temp_path = os.path.join(mh_path, "mh.temp")
        mh_com_path = os.path.join(mh_path, "mh.com")

        if os.path.exists(mh_temp_path):
            os.remove(mh_temp_path)

        if os.path.exists(mh_com_path):
            os.remove(mh_com_path)
         
        high_level_write("mh.temp", [], qR_charges, os.path.join(here, "mh.template"))
        mh = calc.setup_calc("mh", high_level)
        os.chdir(here)
    
        #run model calculations
        #print("Running initial gaussian calculation")
        if low_level == "gaussian":
            ml_proc = ml.run(mol_aug_model, nprocs)
            ml_proc.wait()
            mh_proc = mh.run(mol_aug_model, nprocs)
            mh_proc.wait()
        elif low_level =="xtb":
            print("Running xTB calculation")
            xtb_command = "xtb aug_model.xyz --scc --molden > xtb.out"
            os.chdir("./ml")
            xtb_proc = subprocess.run(xtb_command, shell=True, check=True)
            #xtb_proc.wait()
            os.chdir(here)
            print("Done xTB caclulation")
            print("Running Gaussian calculation")
            mh_proc = mh.run(mol_aug_model, nprocs)
            mh_proc.wait()
            print("Done Gaussian")


        # mh RESP fit
        print("Obtaining formatted checkpoint files for mh calculation:")
        if high_level == "gaussian":
            subprocess.call(['formchk', 'mh/mh.chk'])
            mwfn_run("auto_resp", "mh/mh.fchk", nprocs)
        else:
            raise NotImplementedError(f"model-high RESP fit not implemented in {high_level}") 

        # ml RESP fit
        print("Obtaining formatted checkpoint files for ml calculation:")
        if low_level == "gaussian":
            subprocess.call(['formchk', 'ml/ml.chk'])
            mwfn_run("auto_resp", "ml/ml.fchk", nprocs)
        elif low_level == "xtb":
            subprocess.run(["mv", "ml/molden.input", "ml/ml.molden"])
            mwfn_run("auto_resp", "ml/ml.molden", nprocs)
        
        
        # get new model charges
        mh_resp_charges = [float(line.split()[4]) for line in open("mh.chg") if line.strip()]  
        ml_resp_charges = [float(line.split()[4]) for line in open("ml.chg") if line.strip()]  

        # get link atom charges
        qLinkH_charges = mh_resp_charges[-len(qB_charges):]
        qLinkL_charges = ml_resp_charges[-len(qB_charges):]
        print("{qLinkH}", qLinkH_charges)
        print("{qLinkL}", qLinkL_charges)

        # update qB as per EE(L1) scheme
        averaged_charges = [
            0.5*(qLinkH+qLinkL) for qLinkH,qLinkL in zip(qLinkH_charges,qLinkL_charges) 
        ]
        print("Averaged charges (EE(L1) scheme):", averaged_charges)
        for i, char in enumerate(qB_charges):
            char.q = averaged_charges[i]

        # get constraints for constrained RESP fit
        make_mwfn_constraint(qM_charges, qB_charges)
        print("constraint_file generated.")

        # make mwfn input including constraints
        make_mwfn_config("auto_constrained", auto_constrained)

        # run new RL resp fit with constraints
        subprocess.run(["mv", real_file, "old_"+real_file])
        mwfn_run("auto_constrained", real_file_in,nprocs)
        print("Finished constrained RESP fit")
        
        # update qR charges 
        new_rl_resp_charges = [float(line.split()[4]) for line in open(real_file) if line.strip()]  
        print("New rl RESP charges:",new_rl_resp_charges)
        new_charges = new_rl_resp_charges[-qR_len:]

        qR_charges_old = qR_charges.copy()
        qR_charges_new = []

        for i, qR_charge in enumerate(qR_charges):
            new_qR_charge = qR_charge.copy()  # Create a copy of the qR_charge object
            new_qR_charge.q = new_charges[i]
            qR_charges_new.append(new_qR_charge)

        print("Old charges\n", qR_charges_old)
        print("New charges\n", qR_charges_new)

        #new_qR_charges.append(qR)
        print("Change in charges")
        delta_qs = []
        for old_char, new_char in zip(qR_charges_old, qR_charges_new):
            deltaq = old_char.q-new_char.q
            delta_qs.append(deltaq)
            #print(f"{old_char.elem, old_char.get_pos()} dQ: {deltaq}")
            #print(f"dQ: {deltaq}")
        print("dQ list", delta_qs)
        ave_deltaq = mean(delta_qs)
        max_deltaq = max(delta_qs)
        qR_charges = qR_charges_new
        print("Largest dQ: ", max_deltaq)
        print("Average dQ:", ave_deltaq)
        
        with open(here+"/scc_ee.out", "a") as output_file:
                output_file.write("\n-----------------------------") 
                output_file.write("\nIteration: " + str(iteration))
                output_file.write("\nAverage dQ: "  + str(round(ave_deltaq,5)))
                output_file.write("\nLargest dQ: "  + str(round(max_deltaq,5)))
                
        #check converged
        if  ave_deltaq <= q_thresh:
            scf_active = False
            with open(here+"/scc_ee.out", "a") as output_file:
                output_file.write("\n-----------------------------") 
                output_file.write("\nSCF CONVERGED")
                end_time = datetime.now()
                output_file.write("\nELAPSED TIME: " + str(end_time - start_time))
                output_file.write("\nENDING TIME: " + str(end_time))

            print("SCF converged")
        else:
            print("Not converged, next iteration")

        
if __name__=='__main__':
    main()
