#!/usr/bin/env python
import os,sys
import numpy as np
import glob
import argparse
from shutil import copyfile
from dynamixsampling import Condition

# Functions based (with permission) on the PyRAI2MD package
# Adapted and modified by Federico J. Hernandez 


def getCondition(text):
    """
    Compiles data contained in 'text' into Condition object

    Parameters
    ----------
    text : list<str>
        A list of lines from an initial conditions file, the lines
        are formatted as "symbol X Y Z mass Vx Vy Vz" for each atom
        in the system. Units are Angstrom for coordinates, AMU for
        masses, and Bohr/a.u. of time for velocities

    Return
    ------
    condition : Condition
        A Condition object containing the atomic coordinates, velocities
        and masses of a set of atoms. All necessary data to initialize a
        nonadiabatic dynamics trajectory

    """
    condition = Condition()

    for line in text:
        split_line = line.split()
        for i in range(1,len(split_line)): 
            split_line[i] = float(split_line[i])

        condition.types.append( split_line[0] )
        condition.coordinates.append( [split_line[1], split_line[2], split_line[3] ])
        condition.velocities.append( [split_line[5], split_line[6], split_line[7]] )
        condition.masses.append( split_line[4] )

    condition.coordinates = np.array( condition.coordinates )
    condition.velocities = np.array( condition.velocities )
    condition.masses = np.array( condition.masses )

    return condition


def read_input(in_file):
    """
    Extract initial condition data from a file into a series of Condition objects

    Parameters
    ----------
    in_file : str
        Name of valid input file containing initial conditions data, should
        be output of 'dynamixsampling.py' script
    
    Return
    ------
    conditions : list<Condition>
        List containing one Condition object for each initial condition stored
        in in_file

    """
    n_atoms = 0
    conditions = []

    with open(in_file, 'r') as rf:
        for line in rf:
            # First read number of atoms from file header
            if "NATOM" in line:
                n_atoms = int( line.split()[1] )
            # After reading header parse conditions
            elif n_atoms > 0 and "Condition" in line:
                condition_text = []
                # Put all text for a single condition into a list
                for i in range(n_atoms):
                    condition_text.append(next(rf))
                # Extract data from text with getCondition, and add the
                # returned Condition object to the return list
                conditions.append( getCondition(condition_text) )

    return conditions


def write_data_to_file(data, file_name, directory=None, header=None, types=None, start_line=0, end_line=None):
    """
    Writes a given set of data line-by-line into a file. If 'directory' is
    provided, the file will be created in that directory. Otherwise it will
    be created in the current working directory

    Parameters
    ----------
    data : list<float>
        A list of data to be written to a file. Each element of the list will
        be written to a new line of the file
    file_name : str
        Name of the file to write data to. Optional path to file can be specified
        with the 'directory' option. File will be overwritten if it exists
    header : list<str>

    types : list<str>

    directory : (optional) str
        Absolute path of directory where the file should be written. 
        Default is current working directory (os.getcwd())

    """
    if end_line is None:
        end_line = len(data)
    elif end_line > len(data):
        sys.exit('natoms_flex + natoms_qm exceeds the amunt of atoms present in the initial conditions file')
 
    directory = directory or '.'
    file_path = os.path.join(directory, file_name)

    with open(file_path, 'w+') as wf:
        if header is not None:
            for line in header:
                wf.write(line + "\n")

        for i in range(start_line, min(end_line, len(data))):
            line = data[i]
            line_string = "{:>12.8f} {:>12.8f} {:>12.8f}".format(line[0], line[1], line[2])
            if types is not None:
                line_string = "{:>2s} ".format(types[i]) + line_string
            wf.write(line_string + "\n")

def make_directory_structure(directory=None,phase=1):
    """
    """
    if directory is None:
        directory = os.getcwd()
    elif not os.path.isdir( directory ):
        os.mkdir(directory)

    if phase == 0:
        list_dirs = ["mh"]
    elif phase == 1:
        list_dirs = ["mh", "ml", "rl"]
       
    for subdir in list_dirs:
        os.mkdir( os.path.join( directory, subdir ) )

def copy_files( directory, phase ):
    """
    """
    if phase == 0:
        list_dirs = ["mh"]
    elif phase == 1:
        list_dirs = ["mh", "ml", "rl"]

    for subdir in list_dirs:
    
        dest_path = os.path.join( directory, subdir )
        src_path = os.path.join( os.getcwd(), subdir )
        for f in os.listdir( src_path ):
            src_file = os.path.join(src_path, f)
            if os.path.isfile( src_file ):
                copyfile( src_file, os.path.join(dest_path, f) )

    if phase == 0:
        fro_files = ["fromage.in"]
    elif phase == 1:
        fro_files = ["fromage.in", "shell*"]


    for fro_file in fro_files:
        if '*' in fro_file:
            match_files = glob.glob(os.path.join(os.getcwd(), fro_file))
            for src_file in match_files:
                filename = os.path.basename(src_file)
                dest_file = os.path.join(directory, filename)
                copyfile(src_file, dest_file)
        else:
            # Handle non-pattern filenames as before
            src_file = os.path.join(os.getcwd(), fro_file)
            dest_file = os.path.join(directory, fro_file)
            copyfile(src_file, dest_file)

def setup_conditions(conditions,phase,flex,natoms_flex,natoms_qm):
    """
    Write data contained in a set of Condition objects to a series of files
    for running dynamics trajectories on each initial condition separately

    Parameters
    ----------
    conditions : list<Condition>
        A list of Condition objects
    phase : int 
        phase = 0 --> Gas phase 
        phase = 1 --> Crystal

    flex : boolean
        True Flexible ONIOM

    natoms_flex : int
        Amount of atoms in the QM' region    

    natoms_qm : int
        Amount of atoms in the QM region
    """
    for i,condition in enumerate(conditions, phase):
        directory = f"TRAJ_{i}"
        directory = os.path.join(os.getcwd(), directory)
        make_directory_structure(directory,phase)
        copy_files( directory, phase )
        if flex:
            header_string = [str(natoms_qm),""]
            write_data_to_file(data = condition.coordinates,
                               file_name = "mol.init.xyz", 
                               directory = directory, 
                               header= header_string, 
                               types = condition.types,
                               end_line = natoms_qm) 
            header_string = [str(natoms_flex),""]
            write_data_to_file(data = condition.coordinates, 
                               file_name = "shell_flex.xyz", 
                               directory = directory, 
                               header = header_string, 
                               types = condition.types,
                               start_line = natoms_qm,
                               end_line = natoms_qm + natoms_flex)
            write_data_to_file(data = condition.velocities,
                               file_name = "velocity", 
                               directory = directory,
                               end_line = natoms_qm + natoms_flex)
        else:
            header_string = [ str(len(condition.coordinates)), "" ]
            write_data_to_file(data = condition.coordinates,
                           file_name = "mol.init.xyz",
                           directory = directory,
                           header = header_string,
                           types = condition.types)
            write_data_to_file(condition.velocities, "velocity", directory)


def main():

    usage="""

    Setup module for NAMD in fromage

    Usage:
      python3 setup_dyn.py -i initconds -p 0 or 1
      0: Gas Phase; 1: Molecular Crystal
      setup_dyn.py -h for help

    """

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Setup script for fromage dynamics')
    parser.add_argument('-i','--ifile', type=str, nargs='?', default='initconds', help='File containing initial conditions, output by dynamixsampling.py')
    parser.add_argument('-p','--phase', type=int, default=1, help="Define the phase of the calculation: 0: gas - 1: crystal (Default = 1)")
    parser.add_argument('-f','--flex', type=bool, default=False, help="Define if using a flexible environment (Default = False)")
    parser.add_argument('-nf','--natoms_flex', type=int, default=0, help="Amount of atoms in the QM' flexible region. (Default = 0)")
    parser.add_argument('-nqm','--natoms_qm', type=int, default=0, help="Amount of atoms in the QM region. This option is only used when flex=True. (Default = 0)")
    # Read name of initconds file from command line arguments
    args = parser.parse_args()
    in_file = args.ifile
    phase = args.phase
    flex = args.flex
    natoms_flex = args.natoms_flex
    natoms_qm = args.natoms_qm

    if not os.path.isfile(in_file):
        sys.exit("Initial conditions file {} does not exist".format(in_file))
    if flex and natoms_qm <= 0:
        sys.exit("A flexible environment is required but the number of qm atoms is 0")

    conditions = read_input(in_file)
    setup_conditions( conditions, phase, flex, natoms_flex, natoms_qm )


if __name__ == '__main__':
    main()
