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
import warnings
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
    if "frozen_at" in settings.keys():
        frozen_at = settings["frozen_at"]
        if frozen_at is not None:
            ip["freeze_atoms"] = [int(num) for num in frozen_at]
        else:
            ip["freeze_atoms"] = None

    return ip

def get_freqs(out_file,eigvals,hess_dim,mu_derivs=None):
    """
    Get frequencies and IR intensities
    """
    # Define units conversion
    au2kg = 9.10939e-31 # from electron mass (au) to kg
    amu2kg = 1.66054e-27 # from amu to kg
    au2J = 4.35975e-18 # from Hartrees to Joules
    au2m = 5.29177e-11 # from bohr radius to m
    c = 2.99792e8 # speed of light in m s-1
    m2cm = 100.

    freq_conv_units = au2J / (amu2kg * au2m**2. * m2cm**2.)
    freq_conv_units /= (4.* np.pi**2. * c**2.)

    freqs = np.zeros_like(eigvals)
    intensities = np.zeros_like(eigvals)

    # Change the basis for the dipole derivatives from Cartesian to the basis of normal modes Q 
    if mu_derivs:
        assert mu_derivs.shape[0] == hess_dim
        mu_derivs_Q = mw_hessian.T.dot(mu_derivs)

    for i in range(eigvals.shape[0]):
        if eigvals[i] < 0.:
            freqs[i] = -1. * np.sqrt( np.abs(eigvals[i]) * freq_conv_units)
            # Write warning
            out_file.write("\n")
            out_file.write("In normal modes analysis, found an imaginary mode\n")
            out_file.write("Mode: {:>14.8f}  Frequency = {:>14.8f} cm-1\n".format(
                i+1, freqs[i]))
        else:
          freqs[i] = np.sqrt( eigvals[i] * freq_conv_units)
        if mu_derivs:
            intensities[i] = np.sum(mu_derivs_Q[i,:]**2.) #* 42.255

    return freqs, intensities

def get_center_of_mass(masses,natoms,coords):
    """
    """
    r_CoM = np.sum(coords * masses.reshape((natoms,1))/np.sum(masses),axis=0)
    return r_CoM

def get_I_tensor(masses,coords,r_CoM):
    """
    Compute the moment of inertia I tensor, diagonalise it and
    return the moments of inertia (diagonal elements) and the
    products of inertia (off diagonal elements)
    """
    I_tensor = np.zeros((3,3))
    for i, j in enumerate(coords - r_CoM):
        I_tensor += masses[i]*(np.sum(j**2.)*np.diag(np.ones(3)) - np.outer(j,j))

    I_mom, I_axis = np.linalg.eigh(I_tensor)

    return I_mom, I_axis

def get_trans_and_rot_modes(ncoords,masses,coords_rot,I_axis):
    """
    Get the translation and rotation modes
    """
    tr_rot = np.zeros((ncoords,6))
    tr_rot[0::3,0] = masses**0.5
    tr_rot[1::3,1] = masses**0.5
    tr_rot[2::3,2] = masses**0.5

    for i, mi in enumerate(masses):
        mij = mi**0.5
        for j in range(3):
            tr_rot[3*i+j,3] = + mij*(coords_rot[i,1]*I_axis[j,2]-coords_rot[i,2]*I_axis[j,1])
            tr_rot[3*i+j,4] = - mij*(coords_rot[i,0]*I_axis[j,2]+coords_rot[i,2]*I_axis[j,0])
            tr_rot[3*i+j,5] = + mij*(coords_rot[i,0]*I_axis[j,1]-coords_rot[i,1]*I_axis[j,0])

    u, s, v = np.linalg.svd(tr_rot, full_matrices=True)

    return u[:, 6:]

def write_dmu(dmu_matrix):
    """
    Write the ONIOM dipole derivatives 
    """

    dmu_dim = dmu_matrix.shape[0]
    with open("gs_oniom_dipole_derivatives.dat","a") as out_file:

        for i in range(dmu_dim):
            dmu_str = "{:12.9f} {:12.9f} {:12.9f}".format(
                dmu_matrix[i,0], dmu_matrix[i,1], dmu_matrix[i,2]) + "\n"
            out_file.write(dmu_str)    

    return

def write_FCclasses(in_name1,in_name2,symbols,coords,ener,grads,hess,nmodes,freqs,masses):
    """
    Write the .fcc input file for a FCCclasses calculation
    """
    file_name = in_name1 #"FCclasses_input.fcc"
    file_name2 = in_name2 #"FCclasses_normal_modes.fcc"
    file_name3 = "masses.dat"
    if os.path.exists(file_name):
        os.remove(file_name)
    natoms = len(symbols)
    with open(file_name,"a") as out_file:
        out_file.write("INFO" + "\n")
        out_file.write(" State file generated from fromage" + "\n")
        out_file.write("\n")
        out_file.write("GEOM      UNITS=ANGS" + "\n")
        out_file.write("     %s" % natoms + "\n")
        out_file.write("Geometry from fromage optimisation in xyz format \n")
        for i in range(len(coords)):
            coord_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                symbols[i], coords[i,0], coords[i,1], coords[i,2]) + "\n"
            out_file.write(coord_str)
        out_file.write("\n")
        out_file.write("ENER      UNITS=AU" + "\n")
        out_file.write("     %s" % float(ener) + "\n")
        out_file.write("\n")
        out_file.write("GRAD      UNITS=AU" + "\n")
        _write_common_format(grads, out_file)
        out_file.write("\n")
        out_file.write("HESS      UNITS=AU" + "\n")
        _write_common_format(hess, out_file)

    with open(file_name2,"a") as nm_file:
        for i in range(len(coords)):
            for j in range(3):
                coord_str = "{:12.8f}".format(coords[i,j])
                nm_file.write(coord_str + "\n")
        for i in range(6,nmodes.shape[1]):
            for j in range(nmodes.shape[0]):
                nmode_str = "{:12.8f}".format(nmodes[j,i])
                nm_file.write(nmode_str + "\n")
        for i in range(6,freqs.shape[0]):
            freq_str = "{:12.8f}".format(freqs[i])
            nm_file.write(freq_str + "\n")

    with open(file_name3, "a") as mass_file:
        for i in range(len(masses)):
            mass_str = "{:12.8f}".format(masses[i])
            mass_file.write(mass_str + "\n")
   
    return

def write_vibrations_MOLDEN(in_name,symbols,coord,freqs,inten,in_modes,rmass,str_mode=0):
    """
     Write the vibrational analysis in the Molden format
    """
    file_name = in_name 
    if os.path.exists(file_name):
        os.remove(file_name)
    #Reshape nmodes from a 3Natomsx3Natoms array to a (Nmodes,natoms,3) array
    natoms = len(symbols)
    nmodes = in_modes.shape[1]
    normal_modes = np.zeros((nmodes,natoms,3))
    freqs = freqs.reshape((-1,1))
    inten = inten.reshape((-1,1))
    rmass = rmass.reshape((-1,1))
    nfreqs = len(freqs) - str_mode
    frequencies = '\n'.join(['%12.6f' % x for x in freqs[str_mode:]])
    intensities = '\n'.join(['%12.6f' % x for x in inten[str_mode:]]) 
    rmasses = '\n'.join(['%12.6f' % x for x in rmass[str_mode:]])

    for i in range(nmodes):
       count = 0
       for j in range(0,in_modes.shape[0],3):
            for k in range(3):
                normal_modes[i,count,k] = in_modes[j+k,i]    
            count += 1

    coords = ''
    for n, c in enumerate(coord):
        x, y, z = c
        symb = symbols[n]
        coords += '%-5s %12.8f %12.8f %12.8f\n' % (symb, x, y, z) 

    vibrations = ''
    for i in range(str_mode,nmodes):
        vibrations += '  vibration      %5s\n%s\n' % (
            i + 1 - str_mode, '\n'.join([' '.join(['%12.8f' % y for y in x]) for x in normal_modes[i]]))

    molden = """ [MOLDEN FORMAT]
  [N_FREQ]
%s
  [FREQ]
%s
  [INT]
%s
  [NATOM]
%s
  [FR-COORD]
%s
  [RMASS]
%s
  [FR-NORM-COORD]
%s

""" % (nfreqs,frequencies,intensities,natoms,coords,rmasses,vibrations)

    with open(file_name,"w") as out_file:
        out_file.write(molden)

    return 


"""
    for i in range(nmodes):
       count = 0
       for j in range(0,in_modes.shape[0],3):
            for k in range(3):
                normal_modes[i,count,k] = in_modes[j+k,i]    
            count += 1

    with open(file_name,"a") as out_file:
        out_file.write(" [MOLDEN FORMAT]" + "\n")
        out_file.write(" [N_FREQ]" + "\n")
        out_file.write(" %s" % int(len(freqs)-6) + "\n")
        out_file.write(" [FREQ]" + "\n")
        for i in range(6,freqs.shape[0]):
            freq_str = "{:11.6f}".format(freqs[i])
            out_file.write(freq_str + "\n")
        out_file.write(" [INT]" + "\n")
        for i in range(6,intensities.shape[0]):
            int_str = "{:10.6f}".format(intensities[i])
            out_file.write(int_str + "\n")
        out_file.write(" [NATOM]" + "\n")
        out_file.write(" %s" % natoms + "\n")
        out_file.write(" [FR-COORD]" + "\n")
        for i in range(len(coords)):
            coord_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                symbols[i], coords[i,0], coords[i,1], coords[i,2]) + "\n"
            out_file.write(coord_str)
        out_file.write(" [RMASS] \n")
        for i in range(6,rmasses.shape[0]):
            rmass_str = "{:11.6f}".format(rmasses[i])
            out_file.write(rmass_str + "\n")
        out_file.write(" [FR-NORM-COORD] \n")
        for mode in range(6,normal_modes.shape[0]):
            out_file.write(" vibration		%s" % (mode + 1 - 6) + "\n")
            for atom in range(normal_modes.shape[1]):
                mode_str = "{:10.6f} {:10.6f} {:10.6f}".format(
                    normal_modes[mode,atom,0], normal_modes[mode,atom,1], normal_modes[mode,atom,2]) + "\n"
                out_file.write(mode_str)
        out_file.close()

    return None

   """

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
    
    mol_atoms = Nmodes.mol_atoms
    flex_atoms = Nmodes.flex_atoms
    fixed_atoms = Nmodes.fixed_atoms
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
    freeze_atoms = Nmodes.freeze_atoms

    all_flex_atoms = mol_atoms + flex_atoms

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

    rl_en_gr = rl.read_out(in_pos,
                           in_mol = all_flex_atoms,
                           in_shell = fixed_atoms,
                           natoms_flex = natoms_flex)
    ml_en_gr = ml.read_out(in_pos[:dim_qm], natoms_flex = natoms_flex)
    if high_level == "gaussian_cas":
        mh_en_gr = mh.read_out(in_pos[:dim_qm],natoms_flex = natoms_flex)[0:3]
    else:
        mh_en_gr = mh.read_out(in_pos[:dim_qm], natoms_flex = natoms_flex)

    rl_hess = rl.read_hessian(in_pos, natoms_flex = natoms_flex)
    ml_hess = ml.read_hessian(in_pos[:dim_qm], natoms_flex = natoms_flex)
    mh_hess = mh.read_hessian(in_pos[:dim_qm], natoms_flex = natoms_flex)

    # combine results
 
    en_combo = rl_en_gr[0] - ml_en_gr[0] + mh_en_gr[0]
    gr_combo = rl_en_gr[1] - ml_en_gr[1] + mh_en_gr[1]
    hess_combo = c_low**2. * (rl_hess - ml_hess) + c_high**2. * mh_hess

    if freeze_atoms is not None:
        for atom in freeze_atoms:
            dim = int(atom*3)
            gr_combo[dim-3:dim] = 0.0
            hess_combo[dim-3:dim,dim-3:dim] = 0.0
            
    # if linker atoms are included, hess_combo has to be defined as:
      # where J is the Jacobian that can be easily defined according to
      # the Morokuma's definition https://doi.org/10.1021/cr5004419    
    # hess_combo = c_low**2 * (rl_hess[1] - Jac.T * ml_en_gr[1] * Jac.T) + c_high**2. *  Jac.T * mh_en_gr[1] x Jac.T
    # hess_combo = c_low**2. * (rl_hess - np.matmul(np.matmul(Jac.T, ml_hess), Jac)) + c_high**2. * np.matmul(np.matmul(Jac.T, mh_hess), Jac)

    #hess_out = hess_combo
        
    # get dipole derivatives

    try:
        rl_dmu = rl.read_mu(in_pos, natoms_flex = natoms_flex)
        ml_dmu = ml.read_mu(in_pos[:dim_qm], natoms_flex = natoms_flex)
    except NotImplementedError:

        out_file.write(f"No dipole derivatives (d_mu) for option {low_level} \n.")
        out_file.write(f"fromage continues without the ONIOM d_mu calculation  \n.")
        rl_dmu = None
        ml_dmu = None

    try:
        mh_dmu = mh.read_mu(in_pos[:dim_qm], natoms_flex = natoms_flex)
    except NotImplementedError:

        out_file.write(f"No dipole derivatives (d_mu) for option {high_level} \n.")
        out_file.write(f"fromage continues without the ONIOM d_mu calculation  \n.")
        mh_dmu = None

    if rl_dmu and ml_dmu and mh_dmu:
        if rl_dmu.all() and ml_dmu.all() and mh_dmu.all():
#    if rl_dmu.any() and ml_dmu.any() and mh_dmu.any():
            dmu = rl_dmu - ml_dmu + mh_dmu
    else:
        dmu = None

    if verbose > 1:
        _write_hessians(out_file,hess_combo,mh_hess,ml_hess,rl_hess)
        if dmu:
#        if dmu.any():
            _write_mu_derivatives(out_file,rl_dmu,ml_dmu,mh_dmu)
    # print some updates in the output

    if read_hessian is None:
        out_file.write("------------------------------\n")
        out_file.write("ONIOM Hessian computed " + "\n")
        end_time = datetime.now()
        out_file.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
        out_file.write("ENDING TIME: " + str(end_time) + "\n")
        out_file.flush()

    return (dmu, en_combo, gr_combo, hess_combo)

def _write_hessians(out_file,hess_out,mh_hess,ml_hess,rl_hess):
    """
    write the components of the ONIOM hessian and the ONIOM hessian
    """
    out_file.write("\n")
    out_file.write("--------------------------------------------\n")
    out_file.write("------Writting the Hessian components-------\n")
    out_file.write("--------------------------------------------\n")
    out_file.write("\n")
    out_file.write("High level model Hessian:\n")
    out_file.write("\n")

    _write_common_format(mh_hess, out_file)
    out_file.write("\n")
    out_file.write("Low level model Hessian:\n")
    out_file.write("\n")

    _write_common_format(ml_hess, out_file)
    out_file.write("\n")
    out_file.write("Low level real Hessian:\n")
    out_file.write("\n")

    _write_common_format(rl_hess, out_file)

    hess_file = open("oniom_hessian.txt", "w", 1)
    _write_common_format(hess_out, hess_file)

    return None

def _write_common_format(array, out_file):
    N = array.shape[0]
#    out_file.write("Dimension of array %s x %s" "\n" % (N,N))
    if len(array.shape) == 2:
#        M = array.shape[1]
        row_elements = []
        for i in range(N):
#            for j in range(M):
            for j in range(i + 1): # Iterate only over the lower triangular part
                form_element = f"{array[i,j]:9.7f} "
                row_elements.append(f"{form_element}")
                if len(row_elements) == 5:
                    out_file.write("   ".join(row_elements) + "\n")
                    row_elements = []
        if row_elements:
            out_file.write(" ".join(row_elements) + "\n")
    else:
        for i in range(0, N, 5):
            slice = array[i:i+5]
            line = ' '.join(f"{num:9.7f}" for num in slice)
            out_file.write(line + '\n')

    return None

def _write_mu_derivatives(out_file,rl_dmu,ml_dmu,mh_dmu):
    """
    """
    out_file.write("\n")
    out_file.write("--------------------------------------------\n")
    out_file.write("-Writting the Dipole Derivatives components-\n")
    out_file.write("--------------------------------------------\n")
    out_file.write("\n")
    out_file.write("High level model Dipole derivatives:\n")
    out_file.write("\n")
    for i in range(mh_dmu.shape[0]):
            mh_dmu_str = "{:10.6f}   {:10.6f}   {:10.6f}".format(
                mh_dmu[i,0], mh_dmu[i,1], mh_dmu[i,2]) + "\n"
            out_file.write(mh_dmu_str)
    out_file.write("\n")
    out_file.write("Low level model Dipole derivatives :\n")
    out_file.write("\n")
    for i in range(ml_dmu.shape[0]):
            ml_dmu_str = "{:10.6f}   {:10.6f}   {:10.6f}".format(
                ml_dmu[i,0], ml_dmu[i,1], ml_dmu[i,2]) + "\n"
            out_file.write(ml_dmu_str)
    out_file.write("\n")
    out_file.write("Low level real Dipole derivatives:\n")
    out_file.write("\n")
    for i in range(rl_dmu.shape[0]):
            rl_dmu_str = "{:10.6f}   {:10.6f}   {:10.6f}".format(
                rl_dmu[i,0], rl_dmu[i,1], rl_dmu[i,2]) + "\n"
            out_file.write(rl_dmu_str)
    
    return None

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
        self.natoms = int(init_dict["natoms"])
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
        at_reparam = init_dict["at_reparam"]
        if at_reparam:
            self.at_reparam = []
            self.at_reparam = [int(x) for x in init_dict["at_reparam"]]
            self.at_reparam = np.array(self.at_reparam)
        else:
            self.at_reparam = None
        if "read_hessian" in init_dict.keys(): 
            self.read_hessian = bool_cast(init_dict["read_hessian"])
        else:
            self.read_hessian = None

        if "freeze_atoms" in init_dict.keys():
            self.freeze_atoms = init_dict["freeze_atoms"]
        else:
            self.freeze_atoms = None

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

        out_file = open(self.out_file,'a+')
         
        ########### Define units conversion ###############
        ang2bohr = 1.8897259886
        ###################################################

        masses = self.M[:,0]
        natoms = len(masses)
        mass_mat = np.repeat(self.M,3)
        mr = np.repeat(self.M,3)

        mu_derivs, ener, grads, hessian = sequence_hess(self)

#        mw_hessian = np.matmul(np.matmul(mass_mat, hessian), mass_mat)
        mw_hessian = hessian / np.outer(mass_mat, mass_mat)**0.5

        #diagonalise the mw_hessian to get all the freqs
        eigvals, mw_nmodes = np.linalg.eigh(mw_hessian) 
        # convert the eigenvalues to frequencies
        rmasses = np.zeros_like(eigvals)

        #get rescaled non mass-weighted normal modes
        resc_nmodes = mw_nmodes / np.outer(mass_mat**0.5, np.ones((mw_nmodes.shape[1],)))

        freqs, IR_int = get_freqs(out_file,eigvals,mw_hessian.shape[0],mu_derivs)

        # Get the reduced mass per mode
        for mode in range(len(eigvals)):
            mode_sq = mw_nmodes[:,mode]**2.
            mode_rshp = mode_sq.reshape(natoms,3).sum(axis=1)
            rmasses[mode] = 1. / np.sum(mode_rshp / masses)

        coords_b = np.array(self.in_pos) * ang2bohr

        Freq_file = "oniom.freq.molden"
        write_vibrations_MOLDEN(Freq_file,self.types,coords_b,freqs,IR_int,resc_nmodes,rmasses,0)

        FC_file1 = "FCclasses_input.fcc"
        FC_file2 = "FCclasses_normal_modes.fcc"
        write_FCclasses(FC_file1,FC_file2,self.types,self.in_pos,ener,grads,
                       mw_hessian,resc_nmodes,freqs,masses[:])

        # Now the traslational and rotational normal modes will be projected out the hessian

        r_CoM = get_center_of_mass(masses,natoms,coords_b)
        
        I_mom, I_axis = get_I_tensor(masses,coords_b,r_CoM)

        # rotates coordinates 
        coords_rot = np.dot((coords_b - r_CoM), I_axis)

        ncoords = int(3*len(coords_b))
        trans_rot = get_trans_and_rot_modes(ncoords,masses,coords_rot,I_axis)

        # Diagonalize the dynamic Hessian
        eigvals_rot, eigvecs_rot = np.linalg.eigh(np.dot(trans_rot.T, np.dot(mw_hessian, trans_rot))) 
        # Project out translations and rotations
        mw_nmodes_rot = np.dot(trans_rot, eigvecs_rot)
        # Rescale normal modes to convert them from mass-wighted to cartesian displacements
        resc_nmodes_rot = mw_nmodes_rot / np.outer(mass_mat**0.5, np.ones((mw_nmodes_rot.shape[1],)))

#        freqs_rot, IR_int = get_freqs(out_file,eigvals_rot,mw_hessian.shape[0],mu_derivs)
        
        # Get the reduced mass per mode
#        for mode in range(len(eigvals_rot)):
#            mode_sq = mw_nmodes_rot[:,mode]**2.
#            mode_rshp = mode_sq.reshape(natoms,3).sum(axis=1)
#            rmasses[mode] = 1. / np.sum(mode_rshp / masses)
#
#        Freq_file = "oniom_rot.freq.molden"
#        write_vibrations_MOLDEN(Freq_file,self.types,coords_b,freqs_rot,IR_int,resc_nmodes_rot,rmasses,0)
#
#        FC_file1 = "FCclasses_input_rot.fcc"
#        FC_file2 = "FCclasses_normal_modes_rot.fcc"
#        write_FCclasses(FC_file1,FC_file2,self.types,self.in_pos,ener,grads,
#                       mw_hessian,resc_nmodes_rot,freqs_rot,masses[:])

        if mu_derivs:
            write_dmu(mu_derivs)
    
        return None
#

    def get_freqs(self,eigvals,hess_dim,mu_derivs=None):
        """
        Get frequencies and IR intensities
        """
        # Define units conversion
        au2kg = 9.10939e-31 # from electron mass (au) to kg
        amu2kg = 1.66054e-27 # from amu to kg
        au2J = 4.35975e-18 # from Hartrees to Joules
        au2m = 5.29177e-11 # from bohr radius to m
        c = 2.99792e8 # speed of light in m s-1
        m2cm = 100.

        freq_conv_units = au2J / (amu2kg * au2m**2. * m2cm**2.)
        freq_conv_units /= (4.* np.pi**2. * c**2.)

        freqs = np.zeros_like(eigvals)
        intensities = np.zeros_like(eigvals)

        # Change the basis for the dipole derivatives from Cartesian to the basis of normal modes Q 
        if mu_derivs:
            assert mu_derivs.shape[0] == hess_dim
            mu_derivs_Q = mw_hessian.T.dot(mu_derivs)
 
        for i in range(eigvals.shape[0]):
            if eigvals[i] < 0.:
                freqs[i] = 1. * np.sqrt( np.abs(eigvals[i]) * freq_conv_units)
                # Write warning
                out_file.write("\n")
                out_file.write("In normal modes analysis, found an imaginary mode\n")
                out_file.write("Mode: {:>14.8f}  Frequency = {:>14.8f} cm-1\n".format(
                    i+1, freqs[i]))
            else:
              freqs[i] = np.sqrt( eigvals[i] * freq_conv_units)
            if mu_derivs:
                intensities[i] = np.sum(mu_derivs_Q[i,:]**2.) #* 42.255
 
        return freqs, intensities


    def get_center_of_mass(self,masses,natoms,coords):
        """
        """
        r_CoM = np.sum(coords * masses.reshape((natoms,1))/np.sum(masses),axis=0)
        return r_CoM
 
    def get_I_tensor(self,masses,coords,r_CoM):
        """
        Compute the moment of inertia I tensor, diagonalise it and
        return the moments of inertia (diagonal elements) and the
        products of inertia (off diagonal elements)
        """
        I_tensor = np.zeros((3,3))
        for i, j in enumerate(coords - r_CoM):
            I_tensor += masses[i]*(np.sum(j**2.)*np.diag(np.ones(3)) - np.outer(j,j))

        I_mom, I_axis = np.linalg.eigh(I_tensor)
    
        return I_mom, I_axis

#        I_tensor[0,0] = np.sum(masses * (coords[:,1]**2. + coords[:,2]**2.))
#        I_tensor[1,1] = np.sum(masses * (coords[:,0]**2. + coords[:,2]**2.))
#        I_tensor[2,2] = np.sum(masses * (coords[:,0]**2. + coords[:,1]**2.))
#        for i in range(2):
#           for j in range(i+1,3):
#               I_tensor[i,j] = -1. * np.sum(masses * (coords[:,i] * coords[:,j]))
#               I_tensor[j,i] = I_tensor[i,j]
#
#        eivals, eivecs = np.linalg.eigh(I_tensor)   
#        return eivals, eivecs

    def get_trans_and_rot_modes(self,ncoords,masses,coords_rot,I_axis):
        """
        Get the translation and rotation modes
        """
        tr_rot = np.zeros((ncoords,6))
        tr_rot[0::3,0] = masses**0.5
        tr_rot[1::3,1] = masses**0.5
        tr_rot[2::3,2] = masses**0.5

        for i, mi in enumerate(masses):
            mij = mi**0.5
            for j in range(3):
                tr_rot[3*i+j,3] = + mij*(coords_rot[i,1]*I_axis[j,2]-coords_xyz[i,2]*I_axis[j,1])
                tr_rot[3*i+j,4] = - mij*(coords_rot[i,0]*I_axis[j,2]+coords_xyz[i,2]*I_axis[j,0])
                tr_rot[3*i+j,5] = + mij*(coords_rot[i,0]*I_axis[j,1]-coords_xyz[i,1]*I_axis[j,0])

        u, s, v = np.linalg.svd(tr_rot, full_matrices=True)

        return u[:, 6:]

