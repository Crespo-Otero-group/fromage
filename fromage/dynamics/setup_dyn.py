#!/usr/bin/env python
import os,sys
import numpy as np
import argparse
from shutil import copyfile
from dynamixsampling import Condition




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


def write_data_to_file(data, file_name, directory=None, header=None, types=None):
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
    if directory is None:
        directory = os.getcwd()

    with open(directory + "/" + file_name, 'w+') as wf:
        if header is not None:
            for line in header:
                wf.write(line + "\n")

        for i,line in enumerate(data):
            line_string = "{:>12.8f} {:>12.8f} {:>12.8f}"
            line_string = line_string.format(line[0], line[1], line[2])
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
        fro_files = ["fromage.in", "shell.xyz"]

    for fro_file in fro_files:
#    for fro_file in ["fromage.in", "shell.xyz"]:
        src_file = os.path.join( os.getcwd(), fro_file )
        dest_file = os.path.join( directory, fro_file )
        copyfile( src_file, dest_file )


def setup_conditions(conditions,phase):
    """
    Write data contained in a set of Condition objects to a series of files
    for running dynamics trajectories on each initial condition separately

    Parameters
    ----------
    conditions : list<Condition>
        A list of Condition objects

    """
    for i,condition in enumerate(conditions, phase):
#        directory = "TRAJ_" + "{:05}".format(i+1)
        directory = f"TRAJ_{i}"
        directory = os.path.join(os.getcwd(), directory)
        make_directory_structure(directory,phase)
        copy_files( directory, phase )
        header_string = [ str(len(condition.coordinates)), "" ]
        write_data_to_file(condition.coordinates, "mol.init.xyz", directory, header_string, condition.types)
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
    parser.add_argument('-p','--phase', type=int, default=1, help="Define the pase of the calculation: 0: gas - 1: crystal (Default = 1)")
    # Read name of initconds file from command line arguments
    args = parser.parse_args()
    in_file = args.ifile
    phase = args.phase

    if not os.path.isfile(in_file):
        sys.exit("Initial conditions file {} does not exist".format(in_file))

    conditions = read_input(in_file)
    setup_conditions( conditions, phase )


if __name__ == '__main__':
    main()
