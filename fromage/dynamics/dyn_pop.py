#!/usr/bin/env python
import os,sys
import numpy as np
import argparse
import colorsys



def parse_arguments():
    """
    Set up argument parser and interpret command line

    """
    # Initialize Argument Parser
    parser = argparse.ArgumentParser(description='Electronic State Populations Post-processing Tool')

    # Add argument options to parser
    parser.add_argument('-n','--noplot', action='store_true', help='Flag to turn off automatic generation of a GNUPlot script')
    parser.add_argument('-o','--output', nargs='?', default='populations.pl', help='Name of GNUPlot script to write')
    parser.add_argument('-i','--input', nargs='?', default='output.chk', help='Shared name of all files to read populations from')

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


def get_files(in_file):
    """
    Find input files in all subdirectories of current directory

    Parameters
    ----------
    in_file : str
        Name of the input file to read

    Returns
    -------
    files : List<str>
        A List of full file paths to each valid input file that exists in
        the subdirectories of the current directory

    """
    parent = os.getcwd()
    files = []

    for d in os.listdir(parent):
        if os.path.isdir(d) and os.path.isfile( os.path.join(d,in_file) ):
            files.append( os.path.join(d, in_file) )

    return files


def average_populations(raw_pop, n_states):
    """
    """
    averaged = []
    averaged_dict = {}
    n_files = len(raw_pop)
    max_time = 0

    for f in raw_pop:
        for step in f:
            time = step[0]
            state = step[1]
            if time in averaged_dict.keys():
                averaged_dict[time][state-1] += 1
            else:
                averaged_dict[time] = np.zeros(n_states)
                averaged_dict[time][state-1] += 1
    
    for step in sorted(list(averaged_dict.keys())):
        step_average = [step]
        total_count = np.sum(averaged_dict[step])
        for state in averaged_dict[step]:
            step_average.append( float(state) / float(total_count) )
        averaged.append(step_average)

    return averaged



def extract_populations(files):
    """
    """
    populations = []
    n_states = 0

    for i,f in enumerate(files):
        file_data = []
        with open(f, 'r') as rf:
            for line in rf:
                if len(line.split()) > 4 and isType(line.split()[5], int) and isType(line.split()[3], float):
                    if i == 0 and n_states == 0:
                        n_states = int((len(line.split()) - 13) / 2)
                        print("NSTATES::: {}".format(n_states))
                    file_data.append([ float(line.split()[3]), int(line.split()[5]) ])
        populations.append(file_data)

    return average_populations(populations, n_states)


def write_output(populations):
    """
    """
    n_states = len(populations[0]) - 1

    with open("pop.out",'w+') as wf:
        wf.write("# {:>16s}".format("Time (fs)"))
        for state in range(n_states):
            wf.write(" {:>16s}".format("Root " + str(state+1)))
        wf.write("\n")

        for step in populations:
            step_string = " "
            for data in step:
                step_string += " {:>16.8f}".format(data)
            step_string += "\n"
            wf.write(step_string)

    return None


def write_plot_script(populations):
    """
    """
    plot_script = """set title "S1"

set xrange [0.000000:{max_x}]
set yrange [0.000000:1.000000]
set xlabel 'Time (fs)'
set ylabel 'Population'

set term svg
set out 'populations.pl.svg'
""".format(max_x = np.max([ x[0] for x in populations ]))
    n_states = len(populations[0]) - 1

    for state in range(n_states):
        if state == 0:
            plot_script += "plot 'pop.out' u 1:2 w l tit \"State 0\" lw 2.5 lc rgbcolor \"{color}\", \\\n".format(color = get_color(state, n_states))
        else:
            plot_script += "     ''        u 1:{n2} w l tit \"State {n}\" lw 2.5 lc rgbcolor \"{color}\", \\\n".format(n = state, n2 = state+2, color = get_color(state, n_states))

    with open("populations.pl",'w+') as wf:
        wf.write(plot_script)
    return None


def main():
    """
    """
    # Parse command line arguments
    args = parse_arguments()
    plot = not args.noplot
    # Generate a list of output files
    file_list = get_files(args.input)
    # Construct time-resolved populations array
    populations = extract_populations(file_list)
    # Write data to output file
    write_output(populations)
    # Write plot script, if desired by user
    if plot:
        write_plot_script(populations)


if __name__ == "__main__":
    main()
