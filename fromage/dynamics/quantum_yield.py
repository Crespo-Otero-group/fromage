#!/usr/bin/env python
import os,sys
import numpy as np
import argparse
import colorsys
import random
from itertools import islice
from tqdm import tqdm

#==============================================================#
#================  Basic Coordinate Functions  ================#
#==============================================================#

def distance(A,B):
    """
    Calculate distance between two sets of Cartesian coordinates

    Parameters
    ----------
    A : list<float>
        Array of X, Y, and Z coordinates of first center
    B : list<float>
        Array of X, Y, and Z coordinates of second center

    Returns
    -------
    dist : float
        Distance between center A and center B

    """
    A = np.array(A)
    B = np.array(B)

    dist = np.sqrt(np.sum(np.square(A-B)))

    return dist


def angle(A,B,C):
    """
    """
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)

    Vab = (A-B)/np.linalg.norm(A-B)
    Vcb = (C-B)/np.linalg.norm(C-B)

    angle = np.arccos(np.dot(Vab,Vcb)) * 180 / np.pi

    return angle


def dihedral(A,B,C,D):
    """
    """
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    D = np.array(D)

    Vab = (A-B)/np.linalg.norm(A-B)
    Vbc = (B-C)/np.linalg.norm(B-C)
    Vcd = (C-D)/np.linalg.norm(C-D)

    N1 = np.cross(Vab,Vbc)
    N2 = np.cross(Vbc,Vcd)

    dotProd = np.dot(N1,N2) / np.linalg.norm(np.dot(N1,N2))

    dihedral = np.arccos(dotProd) * 180 / np.pi

    return dihedral


def planeplane(A,B,C,D,E,F):
    """
    """
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    D = np.array(D)
    E = np.array(E)
    F = np.array(F)

    Vab = (A-B)/np.linalg.norm(A-B)
    Vbc = (B-C)/np.linalg.norm(B-C)
    Vde = (D-E)/np.linalg.norm(D-E)
    Vef = (E-F)/np.linalg.norm(E-F)

    N1 = np.cross(Vab,Vbc)
    N2 = np.cross(Vde,Vef)

    dotProd = np.dot(N1,N2) / np.linalg.norm(np.dot(N1,N2))

    planeplane = np.arccos(dotProd) * 180 / np.pi

    return planeplane

#===================================================================#
#=================== Coordinate Class Definition ===================#
#===================================================================#

class Coordinate:
    def __init__(self,coord_type,atom_list,cutoff):
        self.atoms = atom_list
        self.ctype = self.getType(coord_type)
        self.cutoff = cutoff

    def getType(self,type_string):
        if type_string is "r":
            return distance
        elif type_string is "a":
            return angle
        elif type_string is "d":
            return dihedral
        elif type_string is "p":
            return planeplane
        #elif type_string is "c":
        #    return centroid
        else:
            return None

    def getTypeString(self):
        if self.ctype == distance:
            return "R"
        elif self.ctype == angle:
            return "A"
        elif self.ctype == dihedral:
            return "D"
        elif self.ctype == planeplane:
            return "P"

    def calculate(self, *args): #A, B, C=None, D=None, E=None, F=None):
        #return self.ctype(A,B,C,D,E,F)
        return self.ctype(*args)

#===================================================================#
#================= Additional Function Definitions =================#
#===================================================================#

def parse_arguments():
    """
    Set up argument parser and interpret command line

    """
    # Initialize Argument Parser
    parser = argparse.ArgumentParser(description='Electronic State Populations Post-processing Tool')

    # Add argument options to parser
    parser.add_argument('-o','--output', nargs='?', default='qyield.dat', help='Name of output data file (default qyield.dat)')
    parser.add_argument('-i','--input', nargs='?', default='geom_mol.xyz', help='Shared name of all files to read geometries from (default geom_mol.xyz')
    parser.add_argument('-c','--coord', nargs='?', default='Geo.inp', help='Name of file containging geometric coordinates to compute, one per line (default Geo.inp)')
    parser.add_argument('-b','--bcycles', nargs='?', type=int, help='Number of bootstrapping cycles (default = 0, ie no bootstrapping)')
    parser.add_argument('-s','--sample', nargs='?', type=int, help='Number of trajectories to sample for bootstrapping (for bootstrapping only)')

    # Return parsed command line arguments
    return parser.parse_args()


def isType(value, ntype):
    """
    """
    try:
        ntype(value)
        return True
    except ValueError:
        return False


def get_files(args):
    """
    Find input files in all subdirectories of current directory

    Parameters
    ----------
    args : Object
        Object contains parsed command line arguments with name of input file
        to read stored in the "input" attribute

    Returns
    -------
    files : List<str>
        A List of full file paths to each valid input file that exists in
        the subdirectories of the current directory

    """
    in_file = args.input
    parent = os.getcwd()
    files = []

    for d in os.listdir(parent):
        if os.path.isdir(d) and os.path.isfile( os.path.join(d,in_file) ):
            files.append( os.path.join(d, in_file) )

    return files


def parse_coord_text(text):
    """
    Convert text string from coordinate file into a Coordinate object

    """
    parsed_text = text.split()
    coord_type = parsed_text.pop(0)
    cutoff = float(parsed_text.pop())
    atom_list = [ int(x) for x in parsed_text ]

    return Coordinate(coord_type,atom_list,cutoff)


def get_coordinates(args):
    """
    Read input file for a list of internal coordinates to compute

    """
    coord_file = args.coord
    coords = []

    with open(coord_file,'r') as rf:
        coords = rf.readlines()

    coord_list = [ parse_coord_text(x) for x in coords ]

    return coord_list


def compute_coords_for_geom(geom_text,coords):
    """
    """
    values = []

    for coordinate in coords:
        atom_text = [ geom_text[x-1] for x in coordinate.atoms ]
        cartesians = [ [float(x[1]),float(x[2]),float(x[3])] for x in atom_text ]
        value = coordinate.calculate(*cartesians)
        values.append(value)

    return values


def get_coords_from_file(filename, coords, natoms):
    """
    """
    all_values = []
    with open(filename,'r') as rf:
        lines = rf.readlines()
    i = -1 * natoms
    end_geom = [ str(x).split() for x in lines[i:] ]

    values = compute_coords_for_geom(end_geom,coords)

    return values

        

#        while True:
#            lines = []
#            lines_gen = islice(rf, 0, natoms+2)
#            for i,line in enumerate(lines_gen):
#                if i > 1:
#                    lines.append(str(line).split())
#            if not lines:
#                break
#            values = compute_coords_for_geom(lines,coords)
#            all_values.append(values)
#    
#    return all_values


def compute_coordinates(file_list,coord_list):
    """
    """
    all_coordinates = []
    n_coords = len(coord_list)
    natoms = 0

    with open(file_list[0],'r') as rf:
        natom_text = rf.readline().strip()
        if isType(natom_text,int):
            natoms = int(natom_text)
        else:
            print("Error reading number of atoms from file {}".format(file_list[0]))

    for f in tqdm(file_list):
        file_coordinates = get_coords_from_file(f, coord_list, natoms)
        all_coordinates.append(file_coordinates)

    return all_coordinates


def write_output(values, coords, args):
    """
    """
    out_file = args.output

    prod_count = 0
    react_count = 0

    with open(out_file,'w+') as wf:
        wf.write("Traj ")
        for k,coord in enumerate(values[0]):
            coord_header = coords[k].getTypeString()
            for atom in coords[k].atoms:
                coord_header += str(atom)
            wf.write(" {:^12s}".format(coord_header))
        wf.write("\n")
        for j,traj in enumerate(values):
            wf.write("{:^4d} ".format(j))
            prod = True
            for i,coord in enumerate(traj):
                wf.write(" {:^12.2f}".format(coord))
                if coords[i].cutoff > 0 and coord > coords[i].cutoff:
                    prod = prod and True
                elif coords[i].cutoff < 0 and coord < (-1 * coords[i].cutoff):
                    prod = prod and True
                else:
                    prod = prod and False
            if prod:
                wf.write("      Prod\n")
                prod_count += 1
            else:
                wf.write( "      React\n")
                react_count += 1

    total = prod_count + react_count
    print("\nTotal trajectories = {}".format(total))
    print("Reactant trajectories = {} ({}%)".format(react_count, react_count * 100. / total))
    print("Product trajectories = {} ({}%)".format(prod_count, prod_count * 100. / total))

    return None


def do_bootstrapping(args):
    """
    Estimate error in quantum yield with bootstrapping method
    """
    data_file = args.output
    n_cycles = args.bcycles
    n_samples = args.sample
    data = []
    percentage = []

    random.seed()

    print("\nBeginning bootstrapping...")

    with open(data_file,'r') as rf:
        for i,line in enumerate(rf):
            if i > 0:
                data.append(line.split())

    for cycle in tqdm(range(n_cycles)):
        count = 0
        sample_data = random.sample(data,n_samples)
        for sample in sample_data:
            traj_type = sample[-1]
            if traj_type == "Prod":
                count += 1
        percentage.append(count * 100. / n_samples)

    percentage = np.array(percentage)

    MEAN = np.mean(percentage)
    SD = np.std(percentage)

    print_string = "\nYield = {:.2f} ".format(MEAN) + u"\u00B1" + " {:.2f} %\n".format(SD)
    print(print_string)
                

def main():
    """
    """
    # Parse command line arguments
    args = parse_arguments()
    # Generate a list of output files
    file_list = get_files(args)
    # Generate a list of internal coordinates to compute
    coordinate_list = get_coordinates(args)
    # Construct time-resolved coordinates array
    values = compute_coordinates(file_list,coordinate_list)
    # Write data to output file
    write_output(values,coordinate_list,args)

    if args.bcycles is not None:
        do_bootstrapping(args)
            


if __name__ == "__main__":
    main()
