#!/usr/bin/env python
### Dynamic sampling module for neural network and Molcas TSH version 1.0
### Jingbai Li Dec 15 2019
### minor fix on wigner sampling Dec 19 2019 Jingbai Li
### major fix on boltzmann sampling Dec 20 2019 Jingbai Li
### added a function to convert the conditions obtained from
### Newton-X to the format required by fromage Feb 2023 Federico Hernandez
### Supports Molcas, g16, ORCA and Turbomole Oct 27 2023 Federico Hernandez

import os,sys,random
import numpy as np
from optparse import OptionParser
from numpy import linalg as la
from periodic_table import Element
from read_frequencies import *
from fromage.io.parse_config_file import bool_cast

class Condition():
    
    def __init__(self):
        self.coordinates = []
        self.velocities = []
        self.masses = []
        self.charges = []
        self.types = []


def g16format(data,dshape):
    ## formatting data
    dlist=[]
    for i in data:
        dlist+=[float(x) for x in i.split()]
    dlist=np.array(dlist)
    dlist=dlist.reshape(dshape)

    return dlist


def LoadNewtonX(input):
    ## This function read .init.newtonx file (NewtonX final_output file) and return all data as a dict
    ## This function doesn't do initial condition sampling

    with open('%s.init.newtonx' % (input),'r') as initcond:
        nx=initcond.read().splitlines()

    start=0
    end=0
    ensemble=[]
    for n,i in enumerate(nx):
        if   'Geometry' in i:
            start=n
        elif 'Velocity' in i:
            end=n
            break

    natom=end-start-1

    for n,i in enumerate(nx):
        if 'Initial condition' in i:
            coord=[i.split() for i in nx[n+2:n+2+natom]]
            coord=np.array(coord)
            xyz=coord[:,[0,2,3,4]]
            mass=coord[:,[5,1]]
            veloc=[i.split() for i in nx[n+3+natom:n+3+natom*2]]
            veloc=np.array(veloc)
            initcond=np.concatenate((xyz,veloc),axis=1)
            initcond=np.concatenate((initcond,mass),axis=1)
            ensemble.append(initcond)

    return ensemble


def LoadXYZ(input):
    ## This function read .init.xyz file (Gen-FSSH.py .init file) and return all data as a dict
    ## This function doesn't do initial condition sampling

    with open('%s.init.xyz' % (input),'r') as initcond:
        gen=initcond.read().splitlines()

    start=0
    end=0
    ensemble=[]
    for n,i in enumerate(gen):
        if 'Init' in i:
            natom=int(i.split()[2])
            initcond=[i.split() for i in gen[n+1:n+1+natom]]
            initcond=np.array(initcond)
            ensemble.append(initcond)

    return ensemble


def Gaussian():
    ## This function generates standard normal distribution variates from a random number in uniform distribution
    ## This function is used for Boltzmann sampling
    ## This function returns a coefficient to update structure or velocity

    u1 = random.uniform(0, 1)
    u2 = random.uniform(0, 1)
    z = np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2) # Box-Muller transformation
    return float(z)


def Boltzmann(sample):
    """
    """
    # expand atoic mass over x y z coordinates
    amass=np.array([np.ones(3)*i for i in amass]) 
    #standard deviation of mass-weighted position
    sigma_Q=np.sqrt(temp*k_to_au)/np.sqrt(freqs*mu_to_hartree)
    #standard deviation of mass-weighted velocity
    sigma_P=np.sqrt(temp*k_to_au)
    #generates update coordinates and momenta pairs Q and P
    Q_P=np.array([[Gaussian(),Gaussian()] for i in freqs])

    # first column is Q
    Q=Q_P[:,0].reshape(nfreq,1)
    # project standard normal distribution back to position space
    Q*=sigma_Q
    # generate identity array to expand Q
    Qvib=np.array([np.ones((natom,3))*i for i in Q])
     # sum mass-weighted position over all modes
    Qvib=np.sum(vib*Qvib,axis=0)
    # un-weight position in Bohr
    Qvib/=np.sqrt(amass*ma_to_amu)
    # cartesian coordinates in Angstrom
    newc=(xyz+Qvib)*bohr_to_angstrom

    # second column is P
    P=Q_P[:,1].reshape(nfreq,1)
    # convert velocity from m/s to Bohr/au
    P*=sigma_P
    # generate identity array to expand P
    Pvib=np.array([np.ones((natom,3))*i for i in P])
    # sum mass-weighted velocity over all modes
    Pvib=np.sum(vib*Pvib,axis=0)
    # un-weight velocity in Bohr/au
    velo=Pvib/np.sqrt(amass*ma_to_amu)

    inicond=np.concatenate((newc,velo),axis=1)
    inicond=np.concatenate((atoms.reshape((-1,1)),inicond),axis=1)
    inicond=np.concatenate((inicond,amass),axis=1)
    inicond=np.concatenate((inicond,achrg),axis=1)

    return inicond


def gauss_function(frequency, temp):
    pass


def Laguerre(n, x):
    """
    Calculate a Laguerre polynomial of order n

    Parameters
    ----------
    n : int
        Order of the polynomial to compute
    x : float
        X-value at which to compute the value of the polynomial
    Return
    ------
    L : float
        Computed value of the polynomial at X

    """
    L = 1 # L=1 when m=0

    for m in range(1, n+1):
        r = 1
        for mm in range(1, m+1):
            r *= float(n - mm + 1) / mm**2
        L += (-1)**m * r * x**m

    return L


def wigner_function(frequency, temp):
    """
    Generates random positions and momenta to find updated coefficients

    Parameters
    ----------
    frequency : float

    temp : float

    Return
    ------
    Q : float

    P : float


    """
    ## This function generates random position Q and momenta P to find uptdate coifficents
    ## This function calls Laguerre to calculate the polynomial
    ## This function returns accepted Q and P
    max_pop = 0.9999
    ex = frequency / (0.69503 * temp) #vibrational temperature: ex=h*c*mu/(kb*T), 0.69503 convert cm-1 to K
    pop = 0
    lvl_pop = []
    n = -1
    while True:
        n += 1
        pop += np.exp(-1 * ex * n) * (1 - np.exp(-1 * ex))   
        lvl_pop.append(pop[0]) #Note here pop is a numpy array, thus pop[0] is the float number
        # Here is how I obtained this equation:
        # calculate partion function, fP=np.exp(ex*-0.5) /( 1 - np.exp(ex*-1) )
        # calculate population, pop=np.exp(-1*ex*(n+0.5))/fP
        #print('wignerfunction:%d %f %f %f %f'%(n,ex,np.exp(-1*ex*n)*(1-np.exp(-1*ex)),pop,max_pop))
        if pop >= max_pop:
            break
    while True:
        random_state = random.uniform(0, pop) # random generate a state
        n = -1
        for i in lvl_pop:                  # but population is not uniformly distributed over states
            n += 1
            if random_state <= i:          # find the lowest state that has more population than the random state
                break
        Q = random.uniform(0, 1) * 10.0 - 5.0
        P = random.uniform(0, 1) * 10.0 - 5.0
        rho2 = 2 * (Q**2 + P**2)
        W = (-1)**n * Laguerre(n, rho2) * np.exp(-0.5 * rho2)
        R = random.uniform(0, 1)
        #print('N: %d Q: %f P: %f W: %f R: %f' % (N,Q,P,W,R))
        if W > R and W < 1:
            #print('N: %d Q: %f P: %f Rho^2: %f W: %f R: %f' % (N,Q,P,rho2/2,W,R))

            break
    
    return float(Q),float(P)


def getSample(data, temp, method):
    """
    Generate a single initial condition from a probability distribution

    Parameters
    ----------
    freq_data : Normal_Modes object
        A Normal_Modes object containing data read from a QM calculation
        of the vibrational frequencies. Normal_Modes object is returned by
        the read_file() function from any File_Reader()
    temp : float
        Temperature at which to sample the distribution, units = Kelvin

    Return
    ------
    init_cond : Condition object
        An instance of a Condition object which contains the displaced
        coordinates and initial velocities for a single initial condition
        sampled from the Wigner distribution

    """
    mu_to_hartree    = 4.55633518e-6
    ma_to_aum        = 1822.88852
    bohr_to_angstrom = 0.52917724
    k_to_au          = 3.16681156e-6

    init_cond = Condition()

    if method == 'boltzmann':
        sigma_Q = np.sqrt( (temp * k_to_au) / (data.frequencies * mu_to_hartree) )
        sigma_P = np.sqrt( temp * k_to_au )
        sample_function = gauss_function
    else:
        sigma_Q = 1 / np.sqrt(data.frequencies * mu_to_hartree * ma_to_aum)
        sigma_P = np.sqrt( data.frequencies * mu_to_hartree / ma_to_aum )
        sample_function = wigner_function

    masses = np.array([ np.ones(3) * i for i in data.atomic_mass ])

    Q_P = np.array([ sample_function(i, temp) for i in data.frequencies ] )
    Q = Q_P[:,0].reshape(data.n_frequencies, 1)
    P = Q_P[:,1].reshape(data.n_frequencies, 1)

    Q *= sigma_Q
    P *= sigma_P

    Q_displacement = np.array([ np.ones(( data.n_atoms, 3)) * i for i in Q ])
    P_displacement = np.array([ np.ones(( data.n_atoms, 3)) * i for i in P ])

    Q_displacement = np.sum( Q_displacement * data.displacements, axis=0 )
    P_displacement = np.sum( P_displacement * data.displacements, axis=0 )

    if method == 'boltzmann':
        Q_displacement /= np.sqrt( masses * ma_to_aum )
        P_displacement /= np.sqrt( masses * ma_to_aum )

    init_cond.coordinates = (data.atom_coords + Q_displacement) * bohr_to_angstrom
    init_cond.velocities = P_displacement
    init_cond.masses = data.atomic_mass
    init_cond.charges = data.atomic_charge
    init_cond.types = data.atom_types

    return init_cond


def Sampling(in_file, n_conditions, random_seed, temp, sample_method, file_format):
    """
    Generate a set of initial conditions from a probability distribution

    Parameters
    ----------
    in_file : str
        Name of a file containing vibrational frequencies and displacements
    n_conditions : int
        Number of initial conditions to generate
    random_seed : int
        A random seed integer
    temp : float
        Temperature at which to sample the distribution, units = Kelvin
    sample_method : str
        Name of the distribution to sample from, options are 'boltzmann' or
        'wigner'
    file_format : str
        String describing the format of the in_file, one of 'molden', 'bagel',
        'g16', 'orca'

    Return
    ------
    initial_conditions : list<Condition>
        List containing one Condition object for each initial condition generated

    """
    # Use random seed, if provided
    if random_seed != -1:
        random.seed(random_seed)

    # Dict of available File_Reader classes
    file_readers = {
        'molden'    : Molden_Reader,
        'bagel'     : Bagel_Reader,
        'g16'       : Gauss_Reader,
        'turbomole' : Turbomole_Reader, 
        'orca'      : Orca_Reader }

    # Read frequency file with File_Reader determined by file_format
    # If file_format is not implemented print an error message
    if file_format not in file_readers.keys():
        sys.exit("Input file format {} not supported".format(file_format))
    else:
        frequency_data = file_readers[file_format]().read_file(in_file)

    initial_conditions = []

    # Generate initial conditions using sampling function stored in Sample
    for i in range(n_conditions):
        initial_conditions.append( getSample( frequency_data, temp, sample_method ) )

    return initial_conditions


def Equilibrium(input,nesmb,iseed,temp,dist,format):
    ## This function recieves input information and read equilibrium geometry
    ## This function use Readdata to call different functions toextract vibrational frequency and mode
    ## This function returns equilibrium geometry
    ## Import this function for external usage

    callsample=['molden','bagel','g16', 'turbomole']  ## these format need to run sampling
    skipsample=['newtonx','xyz']         ## these format read sampled initial conditions

    ## read in function dictionary
    Readdata={
    'molden':    ReadMolden,
    'bagel':     ReadBagel,
    'g16':       Readg16,
    'newtonx':   LoadNewtonX,
    'turbomole': ReadTurbomole,
    'xyz':     LoadXYZ
    }

    if format in callsample:
        sample=Readdata[format](input)
        atoms=sample['atoms']
        xyz=sample['xyz']*0.529177249
        amass=sample['amass']
        achrg=sample['achrg']
        eqcond=np.concatenate((xyz,np.zeros(xyz.shape)),axis=1)
        eqcond=np.concatenate((atoms.reshape((-1,1)),eqcond),axis=1)
        eqcond=np.concatenate((eqcond,amass),axis=1)
        eqcond=np.concatenate((eqcond,achrg),axis=1)
        return eqcond
    else:
        print('Nothing read from %s' % (input))

def get_NXconds():
    """
    Read the initial condistions previously obtained with Newton-X 
    This function must be called in the same directory where the 
    NX trajectory files TRAJ are present. Only these files must be 
    present, otherwise the function will give an error
    """
    bohr_to_angstrom = 0.52917724

    initial_conditions = []
#    init_cond = Condition()

    all_dirs = next(os.walk('.'))[1]
    for dire in all_dirs:
        init_cond = Condition()
        types, xyz, charges, mass, vel = get_NX_data(dire)
        init_cond.types = types[:,0]
        init_cond.coordinates = xyz.astype(float) * bohr_to_angstrom
        init_cond.charges = charges.astype(float)
        init_cond.masses = mass.astype(float)
        init_cond.velocities = vel.astype(float)
        initial_conditions.append(init_cond)

    return initial_conditions

def get_NX_data(path):
    """
    """
    with open(path+'/geom') as nx_geom_data:
        nx_geom = nx_geom_data.readlines()
    geom = [line.split() for line in nx_geom[:]]
    geom = np.array(geom)
    types = geom[:,[0]]
    charges = geom[:,[1]]
    xyz = geom[:,[2,3,4]]  
    #xyz = geom[:,[0,2,3,4]]
    mass = geom[:,[5]]
    nx_geom_data.close()
    with open(path+'/veloc') as nx_vel_data:
        nx_vel = nx_vel_data.readlines()
    vel = [line.split() for line in nx_vel[:]]
    vel = np.array(vel)
    nx_vel_data.close()

    return types, xyz, charges, mass, vel

def print_conditions(conditions):
    """
    """
    with open("initconds", 'w') as wf:
        wf.write("Initial Conditions File\n")
        wf.write("NATOM  {}\n".format( len(conditions[0].coordinates) ))
        wf.write("NCOND  {}\n".format( len(conditions) ))


        for i,cond in enumerate(conditions):
            wf.write("Condition {}\n".format(str(i)))
            for j,atom in enumerate(cond.coordinates):
                atom_string = "{symb} {x:>12.8f} {y:>12.8f} {z:>12.8f} {mass:>12.8f} {vx:>12.8f} {vy:>12.8f} {vz:>12.8f} \n"
                atom_string = atom_string.format(symb = cond.types[j], 
                                                 x = atom[0], 
                                                 y = atom[1], 
                                                 z = atom[2], 
                                                 mass = cond.masses[j][0], 
                                                 vx = cond.velocities[j][0], 
                                                 vy = cond.velocities[j][1], 
                                                 vz = cond.velocities[j][2] )
                wf.write(atom_string)


def main():
    ## This is the main function 
    ## This function calls Sampling to test the methods.

    usage="""

    Dynamic sampling module for NAMD in fromage

    Usage:
      python3 dynamixsampling.py -i input.freq.molden
      python3 dynamixsampling.py -x 1
      python3 dynamixsampling.py -h for help

    """
    description='Dynamic sampling module for NAMD in fromage'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option('-i',dest='input',type=str,nargs=1,help='Input freq.molden file name')
    parser.add_option('-n',dest='ncond',type=int,nargs=1,help='Number of initial conditions to generate (default = 1)',default=1)
    parser.add_option('-s',dest='iseed',type=int,nargs=1,help='Random number seed (0 .. +inf) (default = random)',default=-1)
    parser.add_option('-t',dest='temp', type=float,nargs=1,help='Temperature in Kelvin (default = 298.15)',default=298.15)
    parser.add_option('-m',dest='method',type=str, nargs=1,help='Sampling method: wigner or boltzmann (default = wigner)',default='wigner')
    parser.add_option('--nx',dest='read_NXconds',type=str,nargs=1,help='Use initial conditions generated with NX (1 to turn it on)',default='0')
    parser.add_option('--md',dest='read_MDconds',type=str,nargs=1,help='Use initial conditions generated from a Molecular dynamics simulation (1 to turn it on)',default='0')
    parser.add_option('--vel',dest='read_velocities',type=str,nargs=1,help='get the velocities from a specific file (1 to turn it on)',default='0')
    parser.add_option('--vf',dest='vel_file',type=str,nargs=1,help='File containing the velocities',default=None)
    parser.add_option('--step',dest='step_size',type=int,nargs=1,help='read the conditions from a MD simulation every step_size steps (default=1)',default=1)

    (options, args) = parser.parse_args()
    read_NXconds = bool_cast(options.read_NXconds)
    read_MDconds = bool_cast(options.read_MDconds)
    read_velocities = bool_cast(options.read_velocities) 
    if options.input == None and read_NXconds == False and read_MDconds == False:
        print (usage)
        exit()

    # Store command line arguments as input parameters
#    if read_NXconds == False:
    input_file = options.input.split('.')[0]
    file_format = options.input.split('.')[-1]
    random_seed = options.iseed
    temperature = options.temp
    sample_method = options.method.lower()
    number_of_conditions = options.ncond

    # run Sampling() with specified parameters
    if read_NXconds:
        ensemble = get_NXconds() 
    elif read_MDconds:
        vel_file = options.vel_file
        step_size = options.step_size
        ensemble = get_MDconds(number_of_conditions,vel_file,step_size)
    else:
        ensemble = Sampling(input_file, number_of_conditions, random_seed, temperature, sample_method, file_format)
    # Print sampled conditions to file
    
    print_conditions(ensemble)
   
 
if __name__ == '__main__':
    main()
