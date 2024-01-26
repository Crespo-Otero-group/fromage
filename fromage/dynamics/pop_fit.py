#!/usr/bin/env python3
import os
import argparse
import random
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import sem
from tqdm import tqdm



def parse_arguments():
    """
    Set up argument parser and interpret command line
    """
    # Initialize Argument Parser
    parser = argparse.ArgumentParser(description='Kinetic Fitting for Electronic Populations')

    # Add argument options
    parser.add_argument('-f','--infile', nargs='?', default='pop.out', help='Name of populations data file (default = pop.out')
    parser.add_argument('-i','--input', nargs='?', default='output.chk', help='Shared name of all files to read populations from (for bootstrapping only, default = output.chk)')
    parser.add_argument('-b','--bcycles', nargs='?', type=int, default=0, help='Number of bootstrapping cycles (default = 0, i.e. no bootstrapping)')
    parser.add_argument('-s','--sample', nargs='?', type=int, default=None, help='Number of trajectories in each bootstrap cycle (for bootstrapping only, default 20 percent of all trajectories)')
    parser.add_argument('-p','--pops', nargs='?', help='List of initial populations in each electronic state, from lowest energy to highest. Ex. 0,0,1 for trajectories initialized in Root 3')
    parser.add_argument('-t','--timestep', nargs='?', type=float, default=0.5, help='Size of time step in the simulations, in femtoseconds (default = 0.5)')

    return parser.parse_args()


def isType(value, ntype):
    """
    """
    try:
        ntype(value)
        return True
    except ValueError:
        return False


def compute_eigenvectors(NSTATES,RATE_GUESS):
    """
    Compute the eigenvalues and eigenvectors for a set of
        rate constants in matrix form
    """

    rate_matrix = np.zeros((NSTATES,NSTATES))

    for i in range(len(RATE_GUESS)):
        j = len(RATE_GUESS) - i
        rate_matrix[j][j] = -1 * RATE_GUESS[i]
        rate_matrix[j-1][j] = RATE_GUESS[i]

    values,vectors = np.linalg.eig(rate_matrix)

    return values,vectors


def compute_scalars(vectors, init_pops):
    """
    Compute scalar coefficients alpha
    """

    return np.dot(np.linalg.inv(vectors), np.array(init_pops))


def compute_populations(values, vectors, alphas, T, DT):
    """
    For every time step compute the state populations using the 
        current model and the guess values of the rates
    """

    populations = []

    for t in np.arange(0,T,DT):
        step_pop = 0

        for state in range(len(values)):
            step_pop += alphas[state][0] * vectors[:,state] * np.exp(values[state]*t)
        if t == 0:
            populations = np.reshape(step_pop,(1,-1))
        else:
            step_pop = np.reshape(step_pop,(1,-1))
            populations = np.concatenate((populations,step_pop), axis=0)

    return populations


def residuals(real, predicted):
    """
    Compute residuals in predicted value, and return as 
        an array of shape (n,)
    """

    real = np.delete(np.array(real), 0, 1).tolist()

    residuals = real - predicted

    residuals = np.array(residuals)
    residuals = np.reshape(residuals,(-1,))

    return residuals


def fun(RATES, real_data, init_pops, T, DT):
    """
    Compute populations from model, and compute residulas
    """
    eigen_vals,eigen_vecs = compute_eigenvectors(len(init_pops), RATES)

    alphas = compute_scalars(eigen_vecs, init_pops)

    populations = compute_populations(eigen_vals, eigen_vecs, alphas, T, DT)

    return residuals(real_data, populations)


def get_data_from_file(file_name, nstates):
    """
    Read population data from file
    """
    populations = []

    with open(file_name,'r') as rf:
        for line in rf:
            if "#" not in line:
                data = [ float(x) for i,x in enumerate(line.split()) if i <= nstates ]
                populations.append(data)

    return np.array(populations)


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
                    file_data.append([ float(line.split()[3]), int(line.split()[5]) ])
        populations.append(file_data)

    return average_populations(populations, n_states)


def get_populations_from_random_sample(in_file, nstates, nsample):
    """
    
    """
    random_subset = []

    list_of_files = get_files(in_file)

    if nsample is not None:
        SAMPLE_SIZE = nsample
    else:
        SAMPLE_SIZE = int(len(list_of_files) / 5)

    for i in range(SAMPLE_SIZE):
        rint = np.random.randint(0, len(list_of_files))
        random_subset.append( list_of_files.pop(rint) )

    subset_populations = extract_populations(random_subset)

    data = []

    for tstep in subset_populations:
        step_data = [ float(x) for i,x in enumerate(tstep) if i <= nstates ]
        data.append(step_data)

    return data


def initialize_populations(pop_list):
    """
    Take string of populations from user and convert to normalized 
        NumPy array
    """
    populations = []
    pop_text = pop_list.split(',')

    for state in pop_text:
        populations.append(float(state))

    # Normalize populations, just in case
    populations = [ x / max(populations) for x in populations ]
    # Convert populations to expected format
    populations = [ [x] for x in populations ]

    return populations


def initialize_rates(nstates):
    """
    Set up a list of guess rates, one for each transition
        between available states
    """
    rates = []

    for process in range(1,nstates):
        rates.append( 1 / (25. * process))

    return rates


def compute_error(values):
    """
    """
    values = [ 1/x for x in values ]
    values = np.array(values)

    MEAN = np.mean(values, axis=0)
    SD = np.std(values, axis=0)
    SEM = sem(values, axis=0)

    return MEAN, SD, SEM
    

def print_error_statistics(mean, sd, sem):
    """
    """

    print("\n############################")
    print("# Bootstrap Time Constants #")
    print("############################\n")

    for i,rate in enumerate(mean):
        start_index = len(mean) - i
        end_index = start_index - 1
        print_string = "S{} --> S{} :: ".format(start_index, end_index)
        if sd[i] == 0:
            print_string += "No population transfer detected"
        else:
            print_string += "{:.2f} ".format(rate) + u"\u00B1" + " {:.2f} fs".format(sd[i])
        print(print_string)
    print("\n")
    return None


def print_fit_results(fit, INIT_RATE_GUESS):
    """
    """
    print("\n#########################")
    print("# Fitted Time Constants #")
    print("#########################\n")
    for i,constant in enumerate(fit.x):
        start_index = len(fit.x) - i
        end_index = start_index - 1
        print_string = "S{} --> S{} ::".format(start_index, end_index)
        if constant == INIT_RATE_GUESS[i]:
            print_string += " No population transfer detected"
        else:
            print_string += " {:.2f} fs".format(1/constant)
        print(print_string)
    print("\n")



def main():
    # Set up argument parser and interpret command line arguments
    args = parse_arguments()

    INIT_POPS = initialize_populations(args.pops)
    INIT_RATE_GUESS = initialize_rates(len(INIT_POPS))
    INPUT_FILE = args.infile
    NSAMPLE = args.sample

    # Read populations data from file
    data = get_data_from_file(INPUT_FILE,len(INIT_POPS))

    # Do non-linear least-squares fitting of rate constants to population data
    # using the matrix methods of: https://doi.org/10.1021/ed067p375
    fit = least_squares(fun, INIT_RATE_GUESS, args=(data, INIT_POPS, len(data) * args.timestep, args.timestep))

    # Print fitted time constants
    print_fit_results(fit, INIT_RATE_GUESS)
    
    # Do bootstrapping, if requested
    if args.bcycles is not None:
        random.seed()
        bootstrap_constants = []

        print("Starting bootstrapping procedure for error estimation...")
        for j in tqdm(range(args.bcycles)):
            r_data = get_populations_from_random_sample(args.input, len(INIT_POPS), NSAMPLE)

            fit = least_squares(fun, INIT_RATE_GUESS, args=(r_data, INIT_POPS, len(r_data) * args.timestep, args.timestep))
            bootstrap_constants.append(fit.x)

        print("Complete!")

        mean,sd,sem = compute_error(bootstrap_constants)

        print_error_statistics(mean, sd, sem)


if __name__ == "__main__":
    main()


