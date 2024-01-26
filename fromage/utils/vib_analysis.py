#!/usr/bin/env python
"""
Minimises the energy or the gap penalty function of a molecule in a cluster.
Also performs Surface Hopping dynamics at CASSCF level with OpenMolcas

Template files are needed as well as an xyz file for the molecule and another
one for the surrounding molecules. Overall the use of subprocess is ugly as it
is repeated 3 or 4 times but it was found to handle memory better than Pool
when interfacing with Gaussian.

Energy units are Hartree inside the program but are printed in eV. Distances are
kept as Angstrom throughout and are converted from Bohr if necessary before
reaching this module

 Written by Federico J Hernandez 20-01-2023
"""
import numpy as np
import subprocess
import os
from datetime import datetime

from fromage.io import read_file as rf
from fromage.utils import array_operations as ao
from fromage.utils import calc
from fromage.io.parse_config_file import bool_cast
from fromage.dynamics.periodic_table import Element

ang2bohr = 1.88973 

####################################################
############# Input Parameters Parser ##############
####################################################

def initNmodesParams(at_symbols,init_pos,settings):
    """
    Parse input parameters and return a dict of settings

    Parameters
    ----------
    at_symbols : List<str>
        Atomic symbols for all atoms in the flexible region (QM+QM')
    init_pos : List<float>
        Initial positions of atoms in the flexible region as [X1, Y1, Z1, X2, Y2, Z2...]
        coordinates units are Angstroms
    settings : Dict
        Unfiltered dict read from fromage.in file

    Returns
    -------
    ip : Dict of trajectory settings

    """
    ip = {}
    ip["types"] = at_symbols
    ip["in_pos"] = np.array(init_pos)
    ip["M"] = np.array([getMass(x) for x in at_symbols])
    ip["natoms"] = len(ip.get("types"))
    ip["low_level"] = settings["low_level"]
    ip["high_level"] = settings["high_level"]
    ip["out_file"] = settings["out_file"]
    if "temp" in settings.keys():
        ip["temp"] = settings["temp"]
    if "verbose" in settings.keys():
        ip["verbose"] = settings["verbose"]
    if "nprocs" in settings.keys():
        ip["nprocs"] = settings["nprocs"]
    if "at_reparam" in settings.keys():
        ip["at_reparam"] = settings["at_reparam"]
    if "scaling_high" in settings.keys():
        ip["scaling_high"] = float(settings["scaling_high"])
    else:
        ip["scaling_high"] = 1.
    if "scaling_low" in settings.keys():
        ip["scaling_low"] = float(settings["scaling_low"])
    else:
        ip["scaling_low"] = 1. 
    if "read_hessian" in settings.keys():
        ip["read_hessian"] = settings["read_hessian"]

    return ip

def write_vibrations_MOLDEN(symbols,coords,freqs,intensities,in_modes):
    """
     Write the vibrational analysis in the Molden format
    """

    #Reshape nmodes froma a 3Natomsx3Natoms array to a (Nmodes,natoms,3) array
    natoms = len(symbols)
    nmodes = in_modes.shape[1]
    normal_modes = np.zeros((nmodes,natoms,3))
    for i in range(nmodes):
       count = 0
       for j in range(0,in_modes.shape[0],3):
            for k in range(3):
                normal_modes[i,count,k] = in_modes[j+k,i]    
            count += 1

#    freqs /= 1.e12
    
    with open("freq.molden","a") as out_file:
        out_file.write(" [MOLDEN FORMAT]" + "\n")
        out_file.write(" [N_FREQ]" + "\n")
        out_file.write(" %s" % len(freqs) + "\n")
        out_file.write(" [FREQ]" + "\n")
        for i in range(freqs.shape[0]):
            freq_str = "{:10.6f}".format(freqs[i])
            out_file.write(freq_str + "\n")
        out_file.write(" [INT]" + "\n")
        for i in range(intensities.shape[0]):
            int_str = "{:10.6f}".format(intensities[i])
            out_file.write(int_str + "\n")
        out_file.write(" [NATOM]" + "\n")
        out_file.write(" %s" % natoms + "\n")
        out_file.write(" [FR-COORD]" + "\n")
        for i in range(len(coords)):
            coord_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                symbols[i], coords[i,0], coords[i,1], coords[i,2]) + "\n"
            out_file.write(coord_str)
        out_file.write(" [FR-NORM-COORD]")
        for mode in range(normal_modes.shape[0]):
            out_file.write(" vibration		%s" % (mode + 1) + "\n")
            for atom in range(normal_modes.shape[1]):
                mode_str = "{:10.6f} {:10.6f} {:10.6f}".format(
                    normal_modes[mode,atom,0], normal_modes[mode,atom,1], normal_modes[mode,atom,2]) + "\n"
                out_file.write(mode_str)
        out_file.close()

    return None
   
def getMass(symbol):
    """
    Returns mass in AMU of element with a given symbol

    Parameters
    ----------
    symbol : str
        Atomic symbol of element 

    Returns
    -------
    mass : float
        Atomic mass in AMU

    """
    return Element(symbol).getMass() 

def get_xyz_cluster(QM_natoms, natoms_flex, flex_atoms, geom_mol = None):
    """
    Function to get the cluster with all the flexible atoms
    (QM + QM') from a previous ONIOM optimization

    Parameters
    ----------
    QM_natoms : int Number of atoms in the QM region

    natoms_flex : int Number of flexible atoms in the 
                      QM' region
    Returns
    -------
    atoms_array : np.array Array with atoms positions
                           QM + QM'(flexible)

    """
    atoms_array = []
    if geom_mol is not None:
        mol_atoms = rf.read_pos("geom_mol.xyz")
        for atom in mol_atoms:
            atoms_array.append([ Element(atom.at_num).getSymbol(),atom.x,atom.y,atom.z ])
    else:
        mol_atoms = rf.read_pos("mol.init.xyz")
        for atom in mol_atoms:
            atoms_array.append([ Element(atom.at_num).getSymbol(),atom.x,atom.y,atom.z ])
        for atom in flex_atoms:
            atoms_array.append([ Element(atom.at_num).getSymbol(),atom.x,atom.y,atom.z ])
    atomic_symbols = [ x[0] for x in atoms_array ]
    atomic_positions = [ [x[1], x[2], x[3]] for x in atoms_array ]
    return (atomic_symbols, atomic_positions)

def sequence_hess(Nmodes):
    """
    Run QM calculations in parallel and write and return results

    This function is a modified version of the sequence() function in
    fro_run_flex.py which computes the Force constants matrix for the 
    model and real systems to solve the ONIOM equation of the Hessian
    for the cluster containing the QM and QM'(flexible) regions.
:
    Parameters
    ----------
    in_pos : list<float>
        Input coordinates in 1D array formi (QM + QM'(flexible))
    mol_atoms : list<Atom>
        List of Atom objects that acts as a template for high layer
    all_atoms: list<Atom>
        List of Atom objects for each atom in the whole cluster
    low_level : str
        Program to use for low level calculation
    high_level : str
        Program to use for high level calculation
    c_high : float
        Freq scaling factor for the high level method
    c_low : float
        Freq scaling factor for the low level method

    Returns
    -------
    hess_out : list<float>
        Hessian matrix in atomic units
        List of spin-orbit coupling between spin states
    """
    
    in_pos = np.array(Nmodes.in_pos).flatten()
    print("")
    print("in_pos in sequence")
    print(in_pos.shape)
    print(in_pos)
    print("")
    
    mol_atoms = Nmodes.mol_atoms
    all_atoms = Nmodes.all_atoms
    QM_natoms = Nmodes.QM_natoms
    natoms_flex = Nmodes.natoms_flex
    dim_qm = Nmodes.dim_qm
    low_level = Nmodes.low
    high_level = Nmodes.high
    at_reparam = Nmodes.at_reparam
    fixed_atoms_array = Nmodes.fixed_atoms_array
    nprocs = Nmodes.nprocs
    c_low = Nmodes.c_low
    c_high = Nmodes.c_high
    read_hessian = Nmodes.read_hessian
    out_file = open(Nmodes.out_file,'a+') 
    verbose = Nmodes.verbose

    rl = calc.setup_calc("rl", low_level)
    ml = calc.setup_calc("ml", low_level)
    mh = calc.setup_calc("mh", high_level)

    # Run the calculations as subprocesses with a maximum of 2 simultameous ones
    # at the same time. This order is optimised for the mh calculation being
    # the longest

    all_pos = np.concatenate((in_pos, fixed_atoms_array), axis = 0)

    if read_hessian is None:
        out_file.write("------------------------------\n")
        out_file.write("Computing ONIOM Hessian" + "\n")
        start_time = datetime.now()
        out_file.write("STARTING TIME: " + str(start_time) + "\n")

        if low_level == "fomo-ci" or low_level == "mopac" and at_reparam is not None:
            rl_proc = rl.run_freq(ao.array2atom(all_atoms, all_pos),nprocs,at_reparam)
            rl_proc.wait()
        else:
            rl_proc = rl.run_freq(atoms = ao.array2atom(all_atoms, all_pos), nprocs = nprocs)
            rl_proc.wait()

        # Get the charges and use them for the model region
        rl_charges_array = rl.read_charges()

        if high_level == "fomo-ci" or high_level == "mopac" and at_reparam is not None:
            mh_proc = mh.run_freq(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                             ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),
                             nprocs, at_reparam)
        else:
            mh_proc = mh.run_freq(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                             ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),nprocs)
            mh_proc.wait()
        if low_level == "fomo-ci" or low_level == "mopac" and at_reparam is not None:
            ml_proc = ml.run_freq(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                             ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]),
                             nprocs,at_reparam)
            ml_proc.wait()
        else:
            ml_proc = ml.run_freq(ao.array2atom(mol_atoms, in_pos[:dim_qm]),
                             ao.array2atom(all_atoms[QM_natoms:], all_pos[dim_qm:], rl_charges_array[QM_natoms:]), nprocs)
            ml_proc.wait()

    # read results. Each x_en_gr is a tuple (energy,gradients,scf_energy)

    rl_hess = rl.read_hessian(in_pos, natoms_flex = natoms_flex)
    ml_hess = ml.read_hessian(in_pos[:dim_qm], natoms_flex = natoms_flex)
    mh_hess = mh.read_hessian(in_pos[:dim_qm], natoms_flex = natoms_flex)

    # combine results

    hess_combo = c_low**2. * (rl_hess - ml_hess) + c_high**2. * mh_hess

    # if linker atoms are included, hess_combo has to be defined as:
      # where J is the Jacobian that can be easily defined according to
      # the Morokuma's definition https://doi.org/10.1021/cr5004419    
    # hess_combo = c_low**2 * (rl_hess[1] - Jac.T * ml_en_gr[1] * Jac.T) + c_high**2. *  Jac.T * mh_en_gr[1] x Jac.T
    # hess_combo = c_low**2. * (rl_hess - np.matmul(np.matmul(Jac.T, ml_hess), Jac)) + c_high**2. * np.matmul(np.matmul(Jac.T, mh_hess), Jac)
      

    mw_hessian = np.matmul(np.matmul(mass_mat, hessian), mass_mat)

    hess_out = hess_combo

    if verbose > 1:
        np.savetxt('hessian',np.column_stack([hess_out]))
        
    # print some updates in the output

    if read_hessian is None:
        out_file.write("------------------------------\n")
        out_file.write("ONIOM Hessian computed " + "\n")
        end_time = datetime.now()
        out_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
        out_file.write("ENDING TIME: " + str(end_time) + "\n")
        out_file.flush()

    return (hess_out)

####################################################
####### Normal modes Initialization Function #######
####################################################

class NormalModes:
    """
    defines a NormalModes object that contains in its attributes all the data necessary for 
    an ONIOM vibrational analysis.
    """

    def __init__(self,init_dict,natoms_flex,mol_atoms,flex_atoms,fixed_atoms,fixed_atoms_array):
        """
        Initialise a NormalModes object with settings provided by the user and parsed by initNmodesParams().
        Below, there is a table with the possible attributes set by the user via fromage.in file, followed 
        y another table of all the parameters used by NormalModes class methods which are initialised
        internally. Therefore, not set by the user.
    
        Parameters
            ----------
            init_dict : Dict
                Dict created by initNmodesParams() function containing trajectory settings
                read from user input in fromage.in file
            natoms_flex : int
                Number of atoms belonging to the flexible QM' region
            mol_atoms : list of atom objects
                Atoms in the inner (QM) region
            flex_atoms : list of atom objects
                Atoms in the intermediate region (QM' flexible)
            fixed_atoms : list of atom objects
                Atoms in the outer region (QM' fixed)
            fixed_atoms_array : np.array
                Array with the positions of the atoms in the outer region  
        """

        self.types = init_dict["types"]
        self.in_pos = init_dict["in_pos"]
        self.M = np.reshape(init_dict["M"],(-1,1))# * 1822.8895
#        self.V = init_dict["V"]
        self.natoms = int(init_dict["natoms"])
#        self.stepSize = float(init_dict["stepsize"])
        self.low = init_dict["low_level"].lower()
        self.high = init_dict["high_level"].lower()
        self.out_file = init_dict["out_file"]
        if "temp" in init_dict.keys():
            self.temp = float(init_dict["temp"])
        else:
            self.temp = 298.15
        if "verbose" in init_dict.keys():
            self.verbose = int(init_dict["verbose"])
        else:
            self.verbose = 1
        if "nprocs" in init_dict.keys():
            self.nprocs = str(init_dict["nprocs"])
        self.c_high = init_dict["scaling_high"]
        self.c_low = init_dict["scaling_low"]
        if "at_reparam" in init_dict.keys():
            self.at_reparam = []
            self.at_reparam = [int(x) for x in init_dict["at_reparam"]]
            self.at_reparam = np.array(self.at_reparam)
        else:
            self.at_reparam = None
        if "read_hessian" in init_dict.keys(): 
            self.read_hessian = bool_cast(init_dict["read_hessian"])
        else:
            self.read_hessian = None

    # Initialize internal trajectory attributes
        self.mol_atoms = mol_atoms
        self.flex_atoms = flex_atoms
        self.fixed_atoms = fixed_atoms
        self.fixed_atoms_array = fixed_atoms_array
        self.all_flex_atoms = mol_atoms + flex_atoms
        self.all_atoms = mol_atoms + flex_atoms + fixed_atoms
        self.QM_natoms = len(mol_atoms)
        self.natoms_flex = natoms_flex
        self.natoms_all_flex = self.QM_natoms + self.natoms_flex
        self.natoms_fixed = len(fixed_atoms)
        self.dim_flex = int(3*natoms_flex)
        self.dim_qm = int(3*self.QM_natoms)
        self.dim_all_flex = int(3*self.natoms_all_flex)
        self.Hess = np.zeros((self.dim_all_flex,self.dim_all_flex))
        self.Freqs = np.zeros(self.dim_all_flex)
        self.Norm_modes = np.zeros((self.dim_all_flex,self.dim_all_flex))
        
    def compute_nmodes(self):
        """
        Function to compute the normal modes and frequencies from an ONIOM Hessian
        """
         
        ########### Define units conversion ###############
        au2kg = 9.10939e-31 # from electron mass (au) to kg
        amu2kg = 1.66054e-27 # from amu to kg
        au2J = 4.35975e-18 # from Hartrees to Joules
        au2m = 5.29177e-11 # from bohr radius to m
        c = 2.99792e8 # speed of light in m s-1
        m2cm = 100.
        ang2bohr = 1.88973
        ###################################################
        
        freq_conv_units = NA * au2J / (amu2kg * au2m**2. * m2cm**2.)

        masses = 1. / np.sqrt(np.repeat(self.M,3))
        mass_mat = np.diag(masses)

        hessian = sequence_hess(self)

        mw_hessian = np.matmul(np.matmul(mass_mat, hessian), mass_mat)

        #diagonalise the mw_hessian to get all the freqs
        eigvals, nmodes_tmp = np.linalg.eigh(mw_hessian) 
        # convert the eigenvalues to frequencies
        freqs = np.zeros_like(eigvals)
        for i in range(eigvals.shape[0]):
            try:
                freqs[i] = np.sqrt( eigvals[i] * freq_conv_units / (4.* np.pi**2. * c**2.) )
            except KeyError:
                freqs[i] = -1. * np.sqrt( np.abs(eigvals[i]) * freq_conv_units / (4.* np.pi**2. * c**2.) )
        intensities = np.zeros_like(freqs)

        # Now the traslational and rotational normal modes will be projected out the hessian
        
        # translate the origin of the system to the center of mass
        r_CM = get_center_of_mass()
        disp_coords = np.zeros_like(self.in_pos)
        for i in range(r_CM.size):
            disp_coords[:,i] = self.in_pos[:,i] - r_CM[i]

        coords_b = self.in_pos * ang2bohr
        write_vibrations_MOLDEN(self.types,coords_b,freqs,intensities,nmodes_tmp)        
    
        return None
#

    def get_center_of_mass(self):
        """
        """
        masses = np.array(self.M)
        tot_mass = np.sum(masses) 
        # convert coord into a np.array and express them in a. u.
        coords = np.array(self.in_pos) * ang2bohr
        CM_pos = np.zeros(3)

        for i in range(3):
            CM_pos[i] = np.sum(masses[:] * coords[:,i]) / tot_mass

        return CM_pos
 
    def get_I_tensor(masses,coords):
        """
        Compute the moment of inertia I tensor, diagonalise it and
        return the moments of inertia (diagonal elements) and the
        products of inertia (off diagonal elements)
        """
        I_tensor = np.zeros((3,3))
        I_tensor[0,0] = np.sum(masses * (coords[:,1]**2. + coords[:,2]**2.))
        I_tensor[1,1] = np.sum(masses * (coords[:,0]**2. + coords[:,2]**2.))
        I_tensor[2,2] = np.sum(masses * (coords[:,0]**2. + coords[:,1]**2.))
        for i in range(2):
           for j in range(i+1,3):
               I_tensor[i,j] = -1. * np.sum(masses * (coords[:,i] * coords[:,j]))
               I_tensor[j,i] = I_tensor[i,j]

        eivals, eivecs = np.linalg.eigh(I_tensor)   
        return eivals, eivecs

    def gen_transformed_coords(masses,coords):
        """
        Coordinates in the rotating and translating frame
    
    
        """
        pass
        return
