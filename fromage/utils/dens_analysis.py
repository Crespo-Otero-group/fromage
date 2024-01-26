#!/usr/bin/env python
# Module written by Federico Hernandez
#      March 2023

import os,sys
import subprocess
import numpy as np
import argparse
from fromage.utils import array_operations as ao
from fromage.dynamics.periodic_table import Element
from fromage.io import read_file as rf
from fromage.io import edit_file as ef
from fromage.io.parse_config_file import bool_cast


def setup_calc(calc_name, calc_type):
    """
    Return a calculation of the correct subclass

    """
    calc_type = calc_type.lower()
    calc_types = {"turbomole" : Turbo_calc}
    try:
        calc = calc_types[calc_type](calc_name)
    except KeyError:
        print("Unercognised program: " + calc_type)

    return calc

class Calc(object):
    """
    Abstract class for calculation objects

    Attributes
    ----------
    calc_name : str
        Name of the calculation, typically rl, ml, mh or mg

    """
    def __init__(self, calc_name_in=None, in_here=os.getcwd()):
        """Constructor which sets the calculation name"""
        self.calc_name = calc_name_in
        self.here = in_here
 
    def run(self, atoms, step_num, point_charges = None, nprocs = None):

        raise NotImplementedError("Please Implement this method")
 
    def post_run(self,options=None):
        "Options to do in the post run"
      
        raise NotImplementedError("Please Implement this method")


def turbo_redefine(atoms):
    """Update Turbomole mos and run actual"""
    FNULL = open(os.devnull, 'w')
    ef.write_coord(atoms)
    # Update mos
    subprocess.call("rm -f mos", shell=True)
    with open("define_feed", "w") as tmp_def_in:
        # define input for Huckel guess
        tmp_def_in.write("\n\n\neht\n\n\n\n\n\n\n\n*\n\n")
    subprocess.call("define < define_feed", stdout=FNULL, shell=True)
    subprocess.call("rm -f define_feed", shell=True)
    subprocess.call("actual -r", shell=True)
    return

class Turbo_calc(Calc):
    """
    Calculation with ADC2 or CC2 with Turbomole 7.0 and 7.6

    """
    def run(self, atoms, step_num, point_charges = None, nprocs = None):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a SH dynamics,
        the placeholder &GRAD has to be included to specify for how 

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Gaussian
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method

        """
        FNULL = open(os.devnull, 'w')

#        os.environ["step_num"] = step_num
        os.environ["nprocs"] = nprocs
        subprocess.call("export PARA_ARCH=SMP", shell=True)
        subprocess.call("export PARNODES=$nprocs", shell=True)        
        dens_dir = 'densities'
        work_dir = 'STEP%s' %(str(step_num))
        dens_dir_path = os.path.join(self.here,dens_dir)
        work_dir_path = os.path.join(dens_dir_path,work_dir)
        subprocess.call("mkdir -p %s/%s" %(dens_dir,work_dir), shell=True)
        subprocess.call("cp JOB_AD/* %s/%s" %(dens_dir,work_dir), shell=True)
#        turbo_path = os.path.join(self.here,work_dir_path)
        os.chdir(work_dir_path)
        turbo_redefine(atoms)
        subprocess.call("cp control old_control", shell=True)
        proc = subprocess.Popen(
            "dscf > dscf.out && ricc2 > ricc2.out", stdout=FNULL, shell=True)        
        os.chdir(self.here)

        return proc

    def post_run(self,step_num,state):
#        FNULL = open(os.devnull, 'w')

        dens_dir = 'densities'
        work_dir = 'STEP%s' %(str(step_num))
        dens_dir_path = os.path.join(self.here,dens_dir)
        work_dir_path = os.path.join(dens_dir_path,work_dir)
        os.chdir(work_dir_path)
        dens_name = "Stp%s_densS%s" %(str(step_num),str(state))
        dens_state_adc = "adcp2-tm0f-1a-00%s.cao" %(str(state))
        dens_state_adc_total = "adcp2-xsdn-1a-00%s-total.cao" %(str(state))
        dens_state_cc = "cc2-tm0f-1a-00%s.cao" %(str(state))
        dens_state_cc_total = "cc2-xsdn-1a-00%s-total.cao" %(str(state))
        with open('old_control') as tmp_file:
            tmp_content = tmp_file.readlines()
        new_content = open('new_control',"w")
        for line in tmp_content:
            if "S1_dens" in line:
                new_content.write(line.replace("S1_dens", dens_name))
            elif "adcp2-xsdn-1a-001-total.cao" in line:
                new_content.write(line.replace("adcp2-xsdn-1a-001-total.cao", dens_state_adc_total))
            elif "adcp2-tm0f-1a-001.cao" in line:
                new_content.write(line.replace("adcp2-tm0f-1a-001.cao", dens_state_adc))
            elif "cc2-xsdn-1a-001-total.cao" in line:
                new_content.write(line.replace("cc2-xsdn-1a-001-total.cao", dens_state_cc_total))
            elif "cc2-tm0f-1a-001.cao" in line:
                new_content.write(line.replace("cc2-tm0f-1a-001.cao", dens_state_cc))
            else:
                new_content.write(line)
        new_content.close()
        subprocess.call("cp new_control control", shell=True)
        proc = subprocess.call("bash fanal.sh", shell=True)
#        os.environ["serial_calc"] = "export OMP_NUM_THREADS=n"
#        subprocess.call("export OMP_NUM_THREADS=n", shell=True)
#        proc = subprocess.Popen("ricc2 -fanal > dens_result.out", shell=True)
        
#        proc = subprocess.Popen(
#            "ricc2 -fanal > dens_result.out", stdout=FNULL,shell=True)
        os.chdir(self.here)
        return proc

    def clean_dir(self,step_num):
        dens_dir = 'densities'
        all_dens_dir = 'all_dens'
        work_dir = 'STEP%s' %(str(step_num))
        dens_dir_path = os.path.join(self.here,dens_dir)
        all_dens_path = os.path.join(dens_dir_path,all_dens_dir)
        work_dir_path = os.path.join(dens_dir_path,work_dir)
        os.chdir(work_dir_path)
        subprocess.call("mv *.cub %s" %(all_dens_path), shell=True)
#        subprocess.call("rm CC*", shell=True)
        os.chdir(self.here)
        return 

def get_densities(in_file,method,start_geom,step_size,all_states,state_file,nprocs):
    """
    Compute the electronic densities
    and NBO population analysis (Optional)
    """
    dens = setup_calc('dens',method)
    all_atoms = rf.read_xyz(in_file)
    all_steps = len(all_atoms)
    dens_dir = 'densities'
    #read the current states
    if state_file != None:
        curr_state = []
        with open(state_file) as f:
            f_content = f.readlines()
        for i, line in enumerate(f_content):
#            if line.strip():
#                if line.split()[0].isdigit():
            curr_state.append(int(line.split()[-1]))
        curr_state = np.array(curr_state)
    else:
        curr_state = None
        states = []
        if len(all_states) > 1:
            states = [int(x) for x in all_states]
        else:
            states = all_states
    
    for step_num in range(start_geom,all_steps,step_size):
        atoms_array = []
        for atom in all_atoms[step_num]:
            atoms_array.append(atom.x)
            atoms_array.append(atom.y)
            atoms_array.append(atom.z)
        atoms_array = np.array(atoms_array)
            
        dens_calc = dens.run(atoms=ao.array2atom(all_atoms[step_num],atoms_array),step_num=step_num,nprocs=nprocs)
        dens_calc.wait()
   
        if curr_state.any() != None:
            dens_cube = dens.post_run(step_num,curr_state[step_num]-1)
        else:
            for state in states:
                dens_cube = dens.post_run(step_num,state)
        dens.clean_dir(step_num)
    return None

def get_nbo_pop(in_file,method,start_geom,step_size,state_file):
    """
    """
    here = os.getcwd()
    reading=False
    dens_dir = 'densities'
    dens_dir_path = os.path.join(here,dens_dir)
    time_step = []
    with open(state_file) as sf:
        sf_content = sf.readlines()
        for i, line in enumerate(sf_content):
            time_step.append(float(line.split()[2]))
        time_step = np.array(time_step)
#
    all_atoms = rf.read_xyz(in_file)
    all_steps = len(all_atoms)
    out_file = open("pop_analysis.dat","w")
    for step_num in range(start_geom,all_steps,step_size):
#        out_file.write(str(step_num) + "\n")
        out_file.write(str(time_step[step_num]) + "\n")
        work_dir = 'STEP%s' %(str(step_num))
        work_dir_path = os.path.join(dens_dir_path,work_dir)
        os.chdir(work_dir_path)
        with open("dens_out.out") as f:
            f_content = f.readlines()
        tmp_lines = []
        last_NBO_pop = len(f_content) - 1 - \
            f_content[::-1].index("Summary of Natural Population Analysis:\n")
        for line in f_content[last_NBO_pop + 5:]:
            if not line.strip()[0].isdigit():  # if line not number
                break
            else:
                tmp_lines.append(line)
        os.chdir(here)
        for j in range(len(tmp_lines)):
            out_file.write(str(tmp_lines[j]))
#
    return None

def run_theodore(in_file,theo_method,start_geom,step_size,state_file):
    """
    """
    here = os.getcwd()
    reading=False
    dens_dir = 'densities'
    dens_dir_path = os.path.join(here,dens_dir)
    with open(state_file) as sf:
        sf_content = sf.readlines()
        for i, line in enumerate(sf_content):
            time_step.append(float(line.split()[2]))
        time_step = np.array(time_step)
#
    out_file = open("Dyn_time.dat","w")
    for step_num in range(start_geom,all_steps,step_size):
#        out_file.write(str(step_num) + "\n")
        out_file.write(str(time_step[step_num]) + "\n")
        work_dir = 'STEP%s' %(str(step_num))
        work_dir_path = os.path.join(dens_dir_path,work_dir)
        os.chdir(work_dir_path)
        os.environ["theo_method"] = theo_method
        subprocess.call("tm2molden",shell=True)
        subprocess.call("theodore $theo_method -f dens_ana.in > theodore_output",shell=True)
        subprocess.call("mkdir theo_outputs",shell=True)
        if theo_method == 'analyze_tden':
            subprocess.call("mv theodore_output *.mld OmFrag.txt ehFrag.txt tden_summ.txt theo_outputs/",shell=True)
        os.chdir(here)
#
    return None


def main():
   ## This is the main function
 
    prog = """
           dens_analysis.py
           """
    epilog="""
    Density analysis module for NAMD or optimizations in fromage
    Supports: Turbomole ADC2/CC2
 
    Usage:
        python3 dens_analysis.py -f input.xyz -m method -c calc_type -i start_geom -s step_size
        python3 dens_analysis.py -h for help
    """
    description='Density analysis module for NAMD or optimizations in fromage'
    
    parser = argparse.ArgumentParser(prog=prog, description=description,epilog=epilog)
    parser.add_argument('-f','--input', type=str, help='Input file name (name.xyz)')
    parser.add_argument('-sf','--state_file', type=str, 
                      help='File with the step and the corresponding current state in th dynamics (name)',default=None)
    parser.add_argument('-m','--method', 
                      type=str, help='EE method selected for the density analysis')
    parser.add_argument('-c','--calc_type', 
                      type=str, help='magnitude to be computed: (orbitals, density, charges, nbo_pop)',default='density')
    parser.add_argument('-i','--start_geom', 
                      type=int, help='The number of the first geometry to start the calculation',default=0)
    parser.add_argument('-s','--step_size', 
                      type=int, help='Step size between the geometries in the input file',default=1)
    parser.add_argument('-n','--states', 
                      type=int, nargs='*', help='Excited state(s) of interest',default=1)
    parser.add_argument('-p','--procs',
                      type=str, help='Number of proces to run parallel calculations',default=1)   
    parser.add_argument('-th','--theodore',
                      type=str, help='Use Theodore3.0 to analyse the results. Also check the -tm for method selection',default=None)
    parser.add_argument('-tm','--theo_method',
                      type=str, help='Analysis method to run Theodore3.0. It requires -th different from "None"',default='analyze_tden')

    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)
    if args.input == None or args.method == None:
        print(epilog)
        exit()

    in_file = args.input
    method = args.method.lower()
    start_geom = args.start_geom
    step_size = args.step_size
    state_file = args.state_file
    states = args.states
    theodore = args.theodore
    theo_method = args.theo_method
#    states = []
#    states = [int(x) for x in args.states]
    nprocs = str(args.procs)
 
    if args.calc_type == 'density':
        if os.path.exists("densities/"):
            subprocess.call("rm -rf densities/", shell=True)
        subprocess.call("mkdir -p densities/all_dens", shell=True)
        
#        subprocess.call("cp JOB_AD/* densities/", shell=True)

        densities = get_densities(in_file,method,start_geom,step_size,states,state_file,nprocs)
    
#
#   Post procesing methods
#    post_methods = ['nbo_pop','theodore']
#    if args.calc_type in post_methods:
    if args.calc_type == 'nbo_pop':
        if os.path.exists("densities/all_dens/"):
            nbo = get_nbo_pop(in_file,method,start_geom,step_size,state_file)
        else:
            print("Missing densities/all_dens/ directory to read the population analysis")
            sys.exit()
    elif args.theodore != None:
        if os.path.exists("densities/all_dens/"):
            theo = run_theodore(in_file,theo_method,start_geom,step_size,state_file)
        else:
            print("Missing densities/all_dens/ directory to read the population analysis")
            sys.exit()
if __name__ == '__main__':
    main()
    
