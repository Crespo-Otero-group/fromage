#!/usr/bin/env python
import os,sys
import numpy as np
import argparse
import colorsys
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
    def __init__(self,coord_type,atom_list):
        self.atoms = atom_list
        self.ctype = self.getType(coord_type)

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
    parser.add_argument('-n','--noplot', action='store_true', help='Flag to turn off automatic generation of a GNUPlot script')
    parser.add_argument('-o','--output', nargs='?', default='ensemble.dat', help='Name of output data file (Default ensemble.dat)')
    parser.add_argument('-i','--input', nargs='?', default='geom_mol.xyz', help='Shared name of all files to read geometries from (Default geom_mol.xyz')
    parser.add_argument('-c','--coord', nargs='?', default='Geo.inp', help='Name of file containging geometric coordinates to compute, one per line')
    parser.add_argument('-d','--data', nargs='?', default='', help='Name of data file from previous motion.py run')
    parser.add_argument('-t','--tstep', nargs='?', type=float, default=0.5, help='Step size of trajectories in femtoseconds (Default=0.5)')
    parser.add_argument('-p','--print', nargs='?', default='1', help='Plot mode. Either a single index of a coordinate to plot as a function of time, or two comma-separated indices to plot two coordinates against each other. Default = 1')

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

def get_color(index, total):
    """
    Return a HEX code for a color that is evenly distributed across the spectrum

    Parameters
    ----------
    index : int
        The index of the color to select from the spectrum
    total : int
        The total number of colors to distribute across the spectrum

    Returns
    -------
    color : str
        A HEX code for a color

    """
    h = float(index) / float(total)
    s = 0.6
    v = 0.8

    r,g,b = colorsys.hsv_to_rgb(h,s,v)

    r = "{:02x}".format(int(r*255))
    g = "{:02x}".format(int(g*255))
    b = "{:02x}".format(int(b*255))

    return "#"+r+g+b


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
    atom_list = [ int(x) for x in parsed_text ]

    return Coordinate(coord_type,atom_list)


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
        while True:
            lines = []
            lines_gen = islice(rf, 0, natoms+2)
            for i,line in enumerate(lines_gen):
                if i > 1:
                    lines.append(str(line).split())
            if not lines:
                break
            values = compute_coords_for_geom(lines,coords)
            all_values.append(values)
    
    return all_values


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

    max_length = len(max(all_coordinates,key = lambda x: len(x)))
    coordinate_values = [ [ [ 0 for z in range(n_coords)] for y in range(max_length)] for x in all_coordinates ]
    coordinate_matrix = [ [ ["NAN" for z in range(n_coords)] for y in range(max_length)] for x in all_coordinates ]
    for i,traj in enumerate(all_coordinates):
        for j,timestep in enumerate(traj):
            for k,coord in enumerate(timestep):
                coordinate_values[i][j][k] = coord
                coordinate_matrix[i][j][k] = coord

    return coordinate_values,coordinate_matrix


def write_output(coord_vals, args):
    """
    """
    out_file = args.output

    with open(out_file,'w+') as wf:
        wf.write("# {:>16s}".format("Time (fs)"))
        for n in range(len(coord_vals[0][0])):
            wf.write(" {:>16s}".format("Coord " + str(n+1)))
        wf.write("\n")

        for step in range(len(coord_vals[0])):
            wf.write("  {:>16.4f}".format(step * args.tstep))
            for traj in coord_vals:
                for coord in traj[step]:
                    if isinstance(coord,float):
                        wf.write(" {:>16.8f}".format(coord))
                    else:
                        wf.write(" {}".format(coord))
            wf.write("\n")

    return None
                

def write_plot_script(coordinates, args):
    """
    Write a GNUPlot script in the user-indicated style

    """
    data_file = args.output
    plottype = [int(x) for x in args.print.split(',')]
    n_coords = len(coordinates[0][0])
    n_traj = len(coordinates)
    min_vals = np.amin(np.amin(np.array(coordinates), axis=1),axis=0)
    max_vals = np.amax(np.amax(np.array(coordinates), axis=1),axis=0)

    script_text = """set title "Ensemble Motion Plot"

set xrange [{XMIN:.6f}:{XMAX:.6f}]
set yrange [{YMIN:.6f}:{YMAX:.6f}]
set xlabel '{XLABEL}'
set ylabel '{YLABEL}'
set key off

set term svg
set out 'motion.pl.svg'

p for [i={START}:{END}:{STEP}] "{FILE}" """

    if len(plottype) > 1:
        xlabel = "Coordinate " + str(plottype[0])
        ylabel = "Coordinate " + str(plottype[1])
        xmin = min_vals[plottype[0]-1]
        xmax = max_vals[plottype[0]-1]
        ymin = min_vals[plottype[1]-1]
        ymax = max_vals[plottype[1]-1]
        start = plottype[0] + 1
        end = ((n_traj - 1)*n_coords) + start
        step = n_coords
        diff = plottype[1] - plottype[0]
        plot_text = "using (column(i)):(column(i+{DIFF})) w l".format(DIFF=diff)
    else:
        xlabel = "Time (fs)"
        ylabel = "Coordinate " + str(plottype[0])
        xmin = 0.0
        xmax = len(coordinates[0]) * args.tstep
        ymin = min_vals[plottype[0]-1]
        ymax = max_vals[plottype[0]-1]
        start = plottype[0] + 1
        end = ((n_traj - 1)*n_coords) + start
        step = n_coords
        plot_text = "u 1:i w l"

    script_text += plot_text

    with open("motion.pl",'w+') as wf:
        wf.write(script_text.format(XMIN = xmin,XMAX = xmax,YMIN=ymin,YMAX=ymax,XLABEL=xlabel,YLABEL=ylabel,START=start,END=end,STEP=step,FILE=data_file))
    return None


def main():
    """
    """
    # Parse command line arguments
    args = parse_arguments()
    plot = not args.noplot
    # Generate a list of output files
    file_list = get_files(args)
    # Generate a list of internal coordinates to compute
    coordinate_list = get_coordinates(args)
    # Construct time-resolved coordinates array
    coord_val,coord_mat = compute_coordinates(file_list,coordinate_list)
    # Write data to output file
    write_output(coord_mat,args)
    # Write plot script, if desired by user
    if plot:
        write_plot_script(coord_val,args)


if __name__ == "__main__":
    main()
