"""Functions to read files from different pieces of chemistry
softare. Most of these functions return lists of Atom objects.
When it comes to energy and gradient values, the units generally
the same as in the file being read. Keep further unit conversion
outside of this file for clarity.
"""
import numpy as np
import fromage.utils.per_table as pt

from fromage.utils.mol import Mol
from fromage.utils import per_table as per
from fromage.utils.atom import Atom
from fromage.utils.dimer import Dimer
from fromage.utils.volume import CubeGrid


def read_vasp(in_name):
    """
    Read VASP POSCAR-like file.

    The real use of this function is to handle one file which contains both
    coordinates and vectors. The actual VASP program is not used anywhere else.
    Make sure the "lattice constant" scaling is set to 1.0, "selective dynamics"
    is not enabled and the file is in Cartesian coordinates.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    M : 3x3 matrix
        Lattice vectors
    atoms : list of Atom types
        Atoms in the file

    """
    with open(in_name) as vasp_file:
        vasp_content = vasp_file.readlines()

    # lattice vectors

    vec1 = vasp_content[2].split()
    vec2 = vasp_content[3].split()
    vec3 = vasp_content[4].split()

    # matrix from vectors
    M = np.zeros((3, 3))
    M[0] = vec1
    M[1] = vec2
    M[2] = vec3

    # reads names of elements and amounts
    species = vasp_content[5].split()
    amounts_str = vasp_content[6].split()
    amounts = map(int, amounts_str)

    # make Atom objects from file
    atoms = []
    for element in species:

        # position of the first and last atom of one kind
        # in the vasp file
        firstAt = 8 + sum(amounts[:species.index(element)])
        lastAt = 8 + sum(amounts[:species.index(element) + 1])

        for line in vasp_content:
            if vasp_content.index(line) in range(firstAt, lastAt):
                xAtom, yAtom, zAtom = map(float, line.split())
                atoms.append(Atom(element, xAtom, yAtom, zAtom))
    return M, atoms


def read_xyz(in_name):
    """
    Read a .xyz file.

    Works for files containing several configurations e.g. a relaxation
    or a trajectory.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    atom_step = list of lists of Atom objects
        Each element of the list represents a configuration of atoms

    """
    with open(in_name) as xyz_file:
        xyz_content = xyz_file.readlines()

    # main list where each element is a relaxation step
    atom_step = []

    for i, line in enumerate(xyz_content):

        # if the line is the amount of atoms in the system
        if line.strip():
            if line.split()[0].isdigit():

                # list of atom objects inside on relaxation step
                atoms = []

                # from 2 lines after the amount of atoms to the last atom line
                # for the relaxation step
                for line_in_step in xyz_content[i + 2:i + int(line) + 2]:
                    elemAtom = line_in_step.split()[0]
                    xAtom, yAtom, zAtom = map(float, line_in_step.split()[1:])
                    atoms.append(Atom(elemAtom, xAtom, yAtom, zAtom))

                atom_step.append(atoms)

    xyz_file.close()
    return atom_step


def read_pos(in_name):
    """
    Return the last or only set of atomic positions in a file

    Currently only .xyz files as they are the most common. To implement more
    types, extend this function by parsing the extension but always return the
    same. read_pos is to be preferred over read_xyz when only one set of
    coordinates is relevant.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    atom_step : list of Atom objects
        The last or only set of atomic positions in the file

    """
    atoms = read_xyz(in_name)[-1]

    return atoms


def mol_from_file(in_name, bonding='', vectors=np.zeros((3, 3))):
    """
    Return a Mol object from a file

    Parameters
    ----------
    in_name : str
        Name of the file to read
    bonding : str
        A string determining the type of bonding in the molecule. Something like
        'dis0.2' or '-13cov'
    vectors : 3 x 3 np array
        The unit cell vectors if pertinent
    Returns
    -------
    mol : Mol object
        The atomic positions in the file as a Mol object

    """
    mol = Mol(read_pos(in_name))
    mol.vectors = vectors
    mol.set_bonding_str(bonding)

    return mol

def dimer_from_file(in_name, bonding=''):
    """
    Return a Dimer object from a file

    The file should have the dimers stated one after the other. The order of
    atoms does not matter as long as they are separated in two.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    bonding : str
        A string determining the type of bonding in the molecule. Something like
        'dis0.2' or '-13cov'
    Returns
    -------
    dim : Dimer object
        The dimer in the file

    """
    double_mol = mol_from_file(in_name)
    mol_a, mol_b = double_mol.split_in_half()
    dim = Dimer(mol_a, mol_b)

    return dim

def read_cp2k(in_name, pop="ESP"):
    """
    Read the charges and energy in a cp2k output file.

    Uses CP2K 4.1 formats. Choose between Mulliken, Hirshfeld or RESP charges.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    pop : str
        Kind of charge to read: mulliken, esp or hirshfeld
    Returns
    -------
    charges : list of floats
        Each partial charge value in the file
    energy : float
        CP2K calculated energy (Kohn-Sham or otherwise) in Hartree

    """
    with open(in_name) as cp2k_file:
        cp2k_content = cp2k_file.readlines()

    if pop.lower() == "mulliken":
        start_tag = "Mulliken Population Analysis"
        char_pos = 4
        line_test = lambda x: x.split()[0].isdigit()
    if pop.lower() == "esp" or pop.lower() == "resp":
        start_tag = " RESP charges:"
        char_pos = 3
        line_test = lambda x: (x.split()[0] == "RESP" and len(x.split()) == 4)
    if pop.lower() in ("hirshfeld", "hirsh"):
        start_tag = "Hirshfeld Charges"
        char_pos = 5
        line_test = lambda x: x.split()[0].isdigit()

    reading = False
    charges = []

    for line in cp2k_content:
        if line.strip():
            if start_tag in line:
                reading = True
            if reading and line_test(line):
                charges.append(float(line.split()[char_pos]))
            if "Total" in line:
                reading = False
            if "ENERGY|" in line:
                energy = float(line.split()[8])

    cp2k_file.close()
    return charges, energy


def read_points(in_name):
    """
    Read point charges from an in-house Ewald.c output.

    The modified version of Ewald.c is needed. The extension is .pts-cry

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    points : Mol object
        Point charges in the file. They have element "point"

    """
    with open(in_name) as pts_file:
        pts_content = pts_file.readlines()

        # store point charges here
    points = Mol([])

    for line in pts_content:
        xIn, yIn, zIn, qIn = map(float, line.split())
        point = Atom("point", xIn, yIn, zIn, qIn)
        points.append(point)

    return points


def read_g_char(in_name, pop="ESP", debug=False):
    """
    Read charges and energy from a Gaussian log file.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    pop : str, optional
        Kind of charge to read, mulliken or esp
    debug : bool, optional
        Return extra energy information. Turn on with care
    Returns
    -------
    charges : list of floats
        Each partial charge value in the file
    energy : float
        Gaussian calculated energy in Hartree
    char_ener : float
        Self energy of the point charges
    n_char : float
        Nuclei-charge interaction energy

    """
    with open(in_name) as gauss_file:
        content = gauss_file.readlines()

    # find last occurrence of Mulliken charges
    if pop.lower() == "mulliken":
        last_mull = len(content) - 1 - \
            content[::-1].index(" Mulliken charges:\n")
    elif pop.lower() == "esp" or pop.lower() == "resp":
        last_mull = len(content) - 1 - \
            content[::-1].index(" ESP charges:\n")
    charges = []

    for line in content[last_mull + 2:]:
        if line.split()[0].isdigit():
            charges.append(float(line.split()[2]))
        else:
            break
    # find each occurrence of Energy
    for line in content:
        if "SCF Done" in line:
            energy = float(line.split()[4])
        if "Total Energy" in line:
            energy = float(line.split()[4])
        if "EIGENVALUE " in line:
            energy = float(line.split()[3])
        if "Self energy of the charges" in line:
            char_ener = float(line.split()[6])
        if "Nuclei-charges interaction" in line:
            n_char = float(line.split()[3])
    if debug:
        return charges, energy, char_ener, n_char
    else:
        return charges, energy


def read_bader(in_name):
    """
    Read charges from a Bader program output file.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    charges : list of floats
        Charges in file

    """
    with open(in_name) as bader_file:
        content_bader = bader_file.readlines()

    # electron charge per atom
    charges = []
    for line in content_bader:
        if line.split()[0].isdigit():
            charge = float(line.split()[4])
            charges.append(charge)

    return charges


def read_qe(in_name):
    """
    Read the final positions of a QE calculation.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    atoms : list of Atom objects
        Last set of atoms in the file

    """
    with open(in_name) as file_qe:
        content = file_qe.readlines()

    last_pos = 0
    for line in content[::-1]:
        if "ATOMIC_POSITIONS" in line.split():
            last_pos = content[::-1].index(line)
            break

    atoms = []
    for line in content[-last_pos:]:
        if line == "End final coordinates\n":
            break
        elem, xPos, yPos, zPos = line.split()
        atom_2_add = Atom(elem, xPos, yPos, zPos, 0)
        atoms.append(atom_2_add)
    return atoms


def read_gauss(in_name):
    """
    Read atoms in a Gaussian input file.

    The format is quite strict, better modify this function before using it.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    atoms : list of Atom objects
        Last set of atoms in the file

    """
    with open(in_name) as file_gauss:
        content = file_gauss.readlines()

    for line in content:
        if line != "\n":
            if line.split()[0].isdigit():
                last_pos = content.index(line)
                break
    atoms = []
    for line in content[last_pos + 1:]:
        if line == "\n" or not line:
            break
        elem, xPos, yPos, zPos = line.split()
        atom_2_add = Atom(elem, xPos, yPos, zPos, 0)
        atoms.append(atom_2_add)
    return atoms


def mol_from_gauss(in_name, pop="ESP"):
    """
    Read Mol from a Gaussian log file

    Include the requested charges.

    Parameters
    ----------
        in_name : str
            Name of the file to read
        pop : str, optional
            Type of charges requested
    Returns
    -------
        out_mol : Mol object
            The atoms in the Gaussian log file, along with their charge

    """
    atoms = read_g_pos(in_name)
    charges = read_g_char(in_name, pop=pop)[0]

    for i, char in enumerate(charges):
        atoms[i].q = char
    out_mol = Mol(atoms)

    return out_mol


def read_fchk(in_name):
    """
    Read a Gaussian .fchk.

    Returns the total energy, gradients and ground state energy.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy : float
        Gaussian total calculated energy in Hartree
    grad : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 etc. Hartree/Bohr
    scf_energy : float
        Gaussian ground state calculated energy in Hartree

    """
    with open(in_name) as data:
        lines = data.readlines()
    grad = []
    reading = False
    for line in lines:
        if line[0].isalpha():
            reading = False
        if reading == True:
            for num in map(float, line.split()):
                grad.append(num)
        if line.startswith("Cartesian Gradient"):
            reading = True
        if line.startswith("Total Energy"):
            energy = float(line.split()[3])
        if line.startswith("SCF Energy"):
            scf_energy = float(line.split()[3])
    grad = np.array(grad)
    return energy, grad, scf_energy

def read_gauss_dyn(in_name,fchk_file,natom,state,states,mult,singlestate,soc_coupling):
    """
    This function is used to read the Gaussian16 TD-DFT output when the nonadiabatic
     dynamics or the Newton-X option is ON.
     Read .log with extended results. 
     Note 1: NACs are only implemeted between the current state in the dynamics
             and S0. 
     Note 2: Reading SOCs from Gaussian is not implemeted yet in fromage. PySOC has to be
             linked.

     Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    ex_energy : float
        Excited state energy in Hartree if any, otherwise the ground state energy
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree
    """
    au2eV = 27.21138386
    energies = []
    ener_temp = []
    gradients = np.array([])
    grad_tmp = []
    nac        = []
    nacs = np.array([])
    nac_tmp = []
    soc        = []

    reading = False
    in_table = False

    with open(fchk_file) as data_fchk:
        lines_fchk = data_fchk.readlines()
    grad = []
    read_g = False
    read_nac = False
    for line in lines_fchk:
        if line[0].isalpha():
            read_g = False
            read_nac = False
        if read_g == True:
            for num in map(float, line.split()):
                grad.append(num)
        if read_nac == True:
            for num in map(float, line.split()):
                nac.append(num)
        if line.startswith("Cartesian Gradient"):
            read_g = True
        if line.startswith("Nonadiabatic coupling"):
            read_nac = True
        if line.startswith("SCF Energy"):
            gr_energy = float(line.split()[3])
#    data_fchk.close()

    with open(in_name) as data:
        lines = data.readlines()

    for line in lines:
        if "Excited State" in line:  
#        if line.startswith("Excited State"):
            ener_temp.append(float(line.split()[4]))

#    data.close()

    # Pack data
    energies.append(gr_energy)
    for i in range(len(ener_temp)):
        energies.append(gr_energy + ener_temp[i]/au2eV)
    energies = np.array(energies)

    gradall = np.zeros((np.sum(states), natom, 3))
    grad = np.array(grad).reshape(natom,3)
    gradall[state - 1] = grad
    gradients = gradall

    ##### I HAVE TO FIX THIS TO COMPLY WITH THE
    ##### ORDER IN THE TEMPLATE FROM NEWTON-X
    nstates = int(np.sum(states))
    ncoup = int(nstates*(nstates-1)/2)
    nacall = np.zeros((ncoup, natom, 3))
    nac = np.array(nac).reshape(natom,3)
    nacall[state - 1] = nac
    nacs = nacall
    soc = np.array(soc)

    return energies, gradients, gr_energy, nac, soc

def read_config(in_name):
    """
    Read a fromage config file.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    settings : dict
        A dictionary with the keywords and the user inputs as strings

    """
    with open(in_name, "r") as data:
        lines = data.readlines()
    settings = {}
    for line in lines:
        if line.strip():
            if line.strip()[0].isalpha():
                if len(line.split()) == 2:
                    settings[line.split()[0].lower()] = line.split()[1]
                else:
                    settings[line.split()[0].lower()] = line.split()[1:]
    return settings


def read_g_pos(in_name):
    """
    Read positions from a Gaussian log file.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    atoms : list of Atom objects
        Atomic positions at the beginning of the file for a single point .log

    """
    with open(in_name) as gauss_file:
        content = gauss_file.readlines()

    # Input orientation
    for i, line in enumerate(content):
        if 'Input orientation:' in line:
            ori_line = i
            break
    atoms = []
    for line in content[ori_line + 5:]:
        if not line.strip()[0].isdigit():  # if line not number
            break
        line_bits = [float(i) for i in line.split()]
        symbol = per.num_to_elem(line_bits[1])
        atom_to_add = Atom(symbol, line_bits[3], line_bits[4], line_bits[5], 0)
        atoms.append(atom_to_add)
    return atoms


def read_ricc2(in_name):
    """
    Read energies and gradients from a Turbomole ricc2.out file.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy : float
        Excited state energy in Hartree if any, otherwise the ground state energy
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    scf_energy : float
        Ground state energy in Hartree

    """
    with open(in_name) as data:
        lines = data.readlines()

    grad_x = []
    grad_y = []
    grad_z = []
    energy = None

    for line in lines:
        if "Total energy of excited state:" in line:
            energy = float(line.split()[5])
        if "Final" in line:
            scf_energy = float(line.split()[5])
        if line.strip():
            if line[0:2] == "dE":
                nums = [float(i.replace("D", "E")) for i in line.split()[1:]]
                if line.split()[0] == "dE/dx":
                    grad_x.extend(nums)
                if line.split()[0] == "dE/dy":
                    grad_y.extend(nums)
                if line.split()[0] == "dE/dz":
                    grad_z.extend(nums)
    grad = []

    # combine in correct format
    for dx, dy, dz in zip(grad_x, grad_y, grad_z):
        grad.append(dx)
        grad.append(dy)
        grad.append(dz)
    # for ground state
    if not energy:
        energy = scf_energy
    grad = np.array(grad)
    return energy, grad, scf_energy

def read_turbo_dyn(in_name, natom, state, states, mult, singlestate, soc_coupling):
    """
     This function is used to read the Turbomole ADC/CC2 output when the nonadiabatic
     dynamics or the Newton-X option is ON.
     Read ricc2 with extended results. 
     Note 1: NACs are not implemeted yet for ricc2 in Turbomole
     Note 2: Reading SOCs from ricc2 is not implemeted yet in fromage

     Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    ex_energy : float
        Excited state energy in Hartree if any, otherwise the ground state energy
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree
    """

    energies = []
    gradients   = np.array([])
    nac        = []
    soc        = []
    grad_x = []
    grad_y = []
    grad_z = []
    grad_tmp = []
 
    reading = False
    in_table = False

    with open(in_name) as data:
        lines = data.readlines()

    for line in lines:
        if "Final" in line:
            gr_energy = float(line.split()[5])
            energies.append(gr_energy)
        if "Energy:" in line:
            energies.append(float(line.split()[1]))
#        if "Total energy of excited state:" in line:
#            energies.append(float(line.split()[5]))
        if reading: 
#        if line.strip():
            if line[0:2] == "dE":
                in_table = True
                nums = [float(i.replace("D", "E")) for i in line.split()[1:]]
                if line.split()[0] == "dE/dx":
                    grad_x.extend(nums)
                if line.split()[0] == "dE/dy":
                    grad_y.extend(nums)
                if line.split()[0] == "dE/dz":
                    grad_z.extend(nums)
        if in_table and "resulting FORCE" in line:
            reading = False
            in_table = False

            # combine in correct format
            for dx, dy, dz in zip(grad_x, grad_y, grad_z):
                grad_tmp.append(dx)
                grad_tmp.append(dy)
                grad_tmp.append(dz)
                
 #           gradients = np.concatenate((gradients, grad_tmp))

            grad_x = []
            grad_y = []
            grad_z = []
        if "cartesian gradient of the energy" in line:
            reading = True 

    if singlestate == 1:
        gradall = np.zeros((np.sum(states), natom, 3))
        grad_tmp = np.array(grad_tmp).reshape(natom,3)
        gradall[state -1] = grad_tmp
        gradients = gradall
    else:
        gradients = np.array(grad_tmp).reshape(np.sum(states), natom, 3)

#    energies = np.reshape(np.array(energies),(-1,1))
    energies = np.array(energies)
    energies[1:] += gr_energy

    nac = np.array(nac)
    soc = np.array(soc)


    return energies, gradients, gr_energy, nac, soc

def read_tbgrad(in_name):
    """
    Read energy gradients from a Turbomole gradient file

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy : float
        Total SCF energy in Hartree
    grad : numpy array
        Energy gradients of the form x1, y1, z1, x2, ... in Hartree/Bohr

    """
    grad = []
    with open(in_name) as data:
        for line in data:
            sline = line.split()
            if "SCF energy" in line:
                energy = float(sline[6])
            if len(sline) == 3 and sline[0] != "$grad":
                for number in sline:
                    grad.append(float(number.replace("D", "E")))

    grad = np.array(grad)
    return energy, grad

def read_tb_MP2_grout(in_name):
    """
    Read energies and gradients from a Turbomole job.last file
    created in a MP2 optimization with just 1 iteration

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy : float
        Excited state energy in Hartree
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    scf_energy : float
        Ground state energy in Hartree

    """

    with open(in_name) as data:
        lines = data.readlines()

    grad_x = []
    grad_y = []
    grad_z = []
    energy = None

    for line in lines:
        if "Total Energy" in line:
            scf_energy = float(line.split()[3])
            energy = scf_energy
        if line.strip():
            if line[0:2] == "dE":
                nums = [float(i.replace("D", "E")) for i in line.split()[1:]]
                if line.split()[0] == "dE/dx":
                    grad_x.extend(nums)
                if line.split()[0] == "dE/dy":
                    grad_y.extend(nums)
                if line.split()[0] == "dE/dz":
                    grad_z.extend(nums)
    grad = []

    # combine in correct format
    for dx, dy, dz in zip(grad_x, grad_y, grad_z):
        grad.append(dx)
        grad.append(dy)
        grad.append(dz)
    # for ground state
    if not energy:
        energy = scf_energy
    grad = np.array(grad)

#    with open(in_name) as data:
#        lines = data.readlines()

#    reading = False

#    for line in lines:
#        if "ENERGY =" in line:
#            energy = float(line.split()[2])
#            scf_energy = energy
#        if "CARTESIAN GRADIENTS" in line:
#            reading = True
#        if reading:
#            if len(line.split()) == 5:
#                if line.split()[-1].isdigit():
#                    nums = [float(i) for i in line.split()[2:]]
#                    grad = np.concatenate((grad, nums))

    return energy, grad, scf_energy

def read_tb_dyn_tddft(in_name, natom, state, states, mult, singlestate, soc_coupling):
    """
     This function is used to read the Turbomole TDDFT output when the nonadiabatic
     dynamics or the Newton-X option is ON.
     Read TDDFT with extended results. 
     Note 1: Reading SOCs from ricc2 is not implemeted yet in fromage

     Parameters
    ----------
    in_name : str
        name of the file to read
    Returns
    -------
    ex_energy : float
        Excited state energy in Hartree if any, otherwise the ground state energy
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree
    """

    energies = []
    gradients   = np.array([])
    nac = []
    soc = []
    grad_x = []
    grad_y = []
    grad_z = []
    grad_tmp = []

    reading = False
    in_table = False


    with open(in_name) as data:
        lines = data.readlines()

    for line in lines:
        if " Total energy:" in line:
            print("Found TDDFT Energies")
            energies.append(float(line.split()[-1]))
        if reading:
            if line[0:2] == "dE":
                in_table = True
                nums = [float(i.replace("D", "E")) for i in line.split()[1:]]
                if line.split()[0] == "dE/dx":
                    grad_x.extend(nums)
                if line.split()[0] == "dE/dy":
                    grad_y.extend(nums)
                if line.split()[0] == "dE/dz":
                    grad_z.extend(nums)
        if in_table and "resulting FORCE" in line:
            reading = False
            in_table = False

            # combine in correct format
            for dx, dy, dz in zip(grad_x, grad_y, grad_z):
                grad_tmp.append(dx)
                grad_tmp.append(dy)
                grad_tmp.append(dz)

            grad_x = []
            grad_y = []
            grad_z = []
        if "cartesian gradient of the energy" in line:
            reading = True

    if singlestate == 1:
        gradall = np.zeros((np.sum(states), natom, 3))
        grad_tmp = np.array(grad_tmp).reshape(natom,3)
        gradall[state -1] = grad_tmp
        gradients = gradall
    else:
        gradients = np.array(grad_tmp).reshape(np.sum(states), natom, 3)

    energies = np.array(energies)
    gr_energy = energies[0]

    nac = np.array(nac)
    soc = np.array(soc)
    return energies, gradients, gr_energy, nac, soc

def read_tb_grout(in_name):
    """
    Read energies and gradients from a Turbomole grad.out TDDFT file.

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy : float
        Excited state energy in Hartree
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    scf_energy : float
        Ground state energy in Hartree

    """
    with open(in_name) as data:
        lines = data.readlines()

    grad_x = []
    grad_y = []
    grad_z = []
    ground_done = False

    for line in lines:
        if "Total energy:" in line:
            if ground_done:
                energy = float(line.split()[2])
            else:
                scf_energy = float(line.split()[2])
                ground_done = True
        if line.strip():
            if line[0:2] == "dE":
                nums = [float(i.replace("D", "E")) for i in line.split()[1:]]
                if line.split()[0] == "dE/dx":
                    grad_x.extend(nums)
                if line.split()[0] == "dE/dy":
                    grad_y.extend(nums)
                if line.split()[0] == "dE/dz":
                    grad_z.extend(nums)
    grad = []

    # combine in correct format
    for dx, dy, dz in zip(grad_x, grad_y, grad_z):
        grad.append(dx)
        grad.append(dy)
        grad.append(dz)
    grad = np.array(grad)

    return energy, grad, scf_energy

def read_xtb(in_name):
    """
    Read energy gradients from a xTB energy and gradient file

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy : float
        Total SCF energy in Hartree
    grad : numpy array
        Energy gradients of the form x1, y1, z1, x2, ... in Hartree/Bohr

    """
    grad = []
    with open(in_name) as data:
        for line in data:
            sline = line.split()
            if "SCF energy" in line:
                try:
                    energy = float(sline[6])
                except ValueError:
                    energy = float(sline[5].split("=")[-1])
            if len(sline) == 3 and sline[0] != "$grad":
                for number in sline:
                    grad.append(float(number.replace("D", "E")))
    grad = np.array(grad)

    return energy, grad

def read_xtb_charges(in_name):
    """
    Read the charges from a xTB  charges file

    Parameters
    ----------
    in_name : str
        Name of the file to read

    Returns
    -------
    array of atom charges 
    
    """
    charges = []
    with open(in_name) as data:
        for line in data:
            charges.append(line.split()[0])
 
    charges = np.array(charges)
 
    return charges

### Functions used in read_molcas_ext

def S2F(M):
    ## This function convert 1D string (e,x,y,z) list to 2D float array
    ## used for read_molcas_ext

    M = [[float(x) for x in row.split()[1: 4]] for row in M]
    return M

def MolcasCoord(M):
    ## This function convert Molcas coordintes to list
    ## used for read_molcas_ext

    coord = []
    for line in M:
        index, a, x, y, z = line.split()[0:5]
        coord.append([a.split(index)[0], float(x), float(y), float(z)])

    return coord

def read_molcas_ext(in_name, natom, state, states, mult, singlestate, soc_coupling):
    """
    This function is used to read the Molcas output when the nonadiabatic dynamics 
    option is ON.
    Read molcas.log with extended results

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    ex_energy : float
        Excited state energy in Hartree if any, otherwise the ground state energy
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree
    """

    print('Molcas read output')

    with open(in_name) as data:
        log = data.readlines()

    spin       = -1
    coord      = []
    casscf     = []
    gradient   = []
    nac        = []
    soc        = []
    soc_mtx    = []
    soc_state  = 0
    sin_state  = 0
    tri_state  = 0


    for i, line in enumerate(log):
        if   """Cartesian coordinates in Angstrom""" in line:
            coord = log[i + 4: i + 4 + natom]
            coord = MolcasCoord(coord)

        elif """Final state energy(ies)""" in line:
            spin += 1
            if   """::    RASSCF root number""" in log[i+3]:
                shift_line = 3  # normal energy output format
                en_col = -1
            else:
                shift_line = 5  # relativistic energy output format
                en_col = 1
            e = [float(x.split()[en_col]) for x in log[i + shift_line: i + shift_line + states[spin]]]
            casscf += e

        elif """Total CASPT2 energies:""" in line:
            if   """::    CASPT2 Root""" in log[i+1]:
                casscf     = []
                shift_line = 1  # normal energy output format
                en_col = -1
            else:
                shift_line = 3  # relativistic energy output format FJH check this!!!!
                en_col = 1
            e = [float(x.split()[en_col]) for x in log[i + shift_line: i + shift_line + states[spin]]]
            casscf += e

        elif """Total MS-CASPT2 energies:""" in line:
            if   """::    MS-CASPT2 Root""" in log[i+1]:
                casscf     = []
                shift_line = 1  # normal energy output format
                en_col = -1
            else:
                shift_line = 3  # relativistic energy output format
                en_col = 1
            e = [float(x.split()[en_col]) for x in log[i + shift_line: i + shift_line + states[spin]]]
            casscf += e

        elif """Total XMS-CASPT2 energies:""" in line:
            if   """::    XMS-CASPT2 Root""" in log[i+1]:
                casscf     = []
                shift_line = 1  # normal energy output format
                en_col = -1
            else:
                shift_line = 3  # relativistic energy output format
                en_col = 1
            e = [float(x.split()[en_col]) for x in log[i + shift_line: i + shift_line + states[spin]]]
            casscf += e

        elif """Molecular gradients """ in line:
            g = log[i + 8: i + 8 + natom]
            g = S2F(g)
            gradient.append(g)

        elif """CI derivative coupling""" in line:
            n = log[i + 8: i + 8 + natom]
            n = S2F(n)
            nac.append(n)

        elif """Nr of states""" in line:
            soc_state = int(line.split()[-1])

        elif """Root nr:""" in line:
            tri_state = int(line.split()[-1])
            sin_state = soc_state - tri_state

        elif """Spin-orbit section""" in line:
            soc_dim = int(sin_state * mult[0] + tri_state * mult[1])
            soc_urt = int(soc_dim * (soc_dim + 1) / 2)
            soc_sfs = np.zeros([soc_dim, soc_dim])
            soc_mtx = np.zeros([soc_state, soc_state])

            # form soc matrix by spin free eigenstates
            for so_el in log[i+11:i+11+soc_urt]:
                i1, s1, ms1, i2, s2, ms2, real_part, imag_part, absolute = so_el.split()
                i1 = int(i1) - 1
                i2 = int(i2) - 1
                va = float(absolute)
                soc_sfs[i1, i2] = va
                soc_sfs[i2, i1] = va

            # reduce soc matrix into configuration state, effective soc
            for s1 in range(sin_state):
                for s2 in range(tri_state):
                    p2 = sin_state + s2
                    first_col = int(sin_state + s2 * mult[1])
                    final_col = int(sin_state + (s2 + 1) * mult[1])
                    soc_mtx[s1, p2] = np.sum(soc_sfs[s1, first_col: final_col]**2)**0.5
                    soc_mtx[p2, s1] = soc_mtx[s1, p2]

    ## extract soc matrix elements
    if len(soc_coupling) > 0 and len(soc_mtx) > 0:
        for pair in soc_coupling:
            s1, s2 = pair
            socme = float(soc_mtx[s1, s2 - states[0] + sin_state]) ## assume low spin is in front of high spin (maybe generalize later)
            soc.append(socme)

    ## pack data
    energy   = np.array(casscf)

    if singlestate == 1:
        gradall = np.zeros((np.sum(states), natom, 3))
        gradall[state - 1] = np.array(gradient)
        gradient = gradall
    else:
        gradient = np.array(gradient)

    nac = np.array(nac)
    soc = np.array(soc)

    gr_energy = energy[0]
    return energy, gradient, gr_energy, nac, soc

def read_molcas(in_name):
    """
    Read energies and gradients from a Molcas .log file with 2 roots.

    Parameters
    ----------
    in_name : str
        Name of the file to read
#    state : int
#        Selected state for which the energy and gradient will ve saved
    Returns
    -------
    ex_energy : float
        Excited state energy in Hartree if any, otherwise the ground state energy
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree

    """
    with open(in_name) as data:
        lines = data.readlines()

    # Initialize cur_state to None; 
    # it will hold the extracted number 
    # if it finds any
    cur_state = None

    for line in lines:
        stripped_line = line.strip().lower()        
        # Check if the line starts with 'rlxroot' after being stripped and lowered
        if stripped_line.startswith('rlxroot'):
            parts = stripped_line.split('=')
            # Check if there is an element after '=' to handle 'Rlxroot=' with no value
            if len(parts) > 1:
                cur_state = parts[1].strip()

    grad = np.array([])
    ex_energy = None
    gr_energy = None

    reading = False
    in_table = False

    for line in lines:
        if line.strip():
            if "RASSCF root number  1 Total energy:" in line:
                gr_energy = float(line.split()[-1])
            if "CASPT2 Root  1     Total energy:" in line:
                gr_energy = None
                gr_energy = float(line.split()[-1])
            if "MS-CASPT2 Root  1     Total energy:" in line:
                gr_energy = None
                gr_energy = float(line.split()[-1])
            if "XMS-CASPT2 Root  1     Total energy:" in line:
                gr_energy = None
                gr_energy = float(line.split()[-1])
            if "RASSCF energy for state" in line:
                ex_energy = float(line.split()[-1])
            if "::    CASPT2 Root" in line:
                words = line.split()
                if len(words) > 5 and words[3] == str(cur_state):
                    ex_energy = float(words[-1])
            if '::    MS-CASPT2 Root' in line:
                words = line.split()
                if len(words) > 5 and words[3] == str(cur_state):
                    ex_energy = float(words[-1])
            if '::    XMS-CASPT2 Root' in line:
                words = line.split()
                if len(words) > 5 and words[3] == str(cur_state):
                    ex_energy = float(words[-1])
            # Gradients
            if "Molecular gradients" in line:
                reading = True
            if reading:
                if len(line.split()) == 4:
                    if len(line.split()[0]) > 1:
                        if line.split()[0][-1].isdigit():
                            nums = [float(i) for i in line.split()[1:]]
                            grad = np.concatenate((grad, nums))
    if not ex_energy:
        ex_energy = gr_energy
    return ex_energy, grad, gr_energy

def read_dftb_out(in_name):
    """
    Read a dftb+ gradient detailed.out file

    Parameters
    ----------
    in_name : str
        name of the file to read
    Returns
    -------
    ex_energy : float
        Excited state energy in Hartree
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree

    """

    grad = []
    exci = 0
    with open(in_name) as lines:
        read_grad = False
        for line in lines:
            if not line.strip():
                read_grad = False
            if read_grad:
#                atom_grads = [float(i) for i in line.split()]
                atom_grads = [float(i) for i in line.split()[1:]] # FJH
                grad.extend(atom_grads)
            if "Total Forces" in line:
                read_grad = True
            if "Total energy" in line:
                gr_energy = float(line.split()[2])
            if "Excitation Energy" in line:
                exci += float(line.split()[2])

    # the detailed out prints forces so *-1 for forces
    grad = -np.array(grad)
    ex_energy = gr_energy + exci

    return ex_energy, grad, gr_energy

def read_dftb_charges(in_name):
    """
    Read a dftb+  Mulliken charges detailed.out file

    Parameters
    ----------
    in_name : str
        name of the file to read

    Returns
    -------
    charges : numpy array of floats
        Charges in the form q1,q2,q3,q4 etc. in Hartree/Bohr

    """

    charges = []
    charges_cm5 = []
    with open(in_name) as lines:
        read_charges = False
        read_charges_cm5 = False
        for line in lines:
            if not line.strip():
                read_charges = False
                read_charges_cm5 = False
            if read_charges:
                charges.append(line.split()[-1])
            # Info about CM5 corrected charges can be
            # found here: J. Chem. Theory Comput. 2012, 8, 2, 527-541
            if read_charges_cm5:
                charges_cm5.append(line.split()[-1])
            if "Atomic gross charges" in line:
                read_charges = True
            if "CM5 corrected atomic gross charges" in line:
                read_charges_cm5 = True

    if len(charges_cm5) == len(charges):
        charges = np.array(charges_cm5)
    else:
        charges = np.array(charges)

    return charges

def read_g_cas(in_name):
    """
    Read a Gaussian .log file for CAS calculations

    Returns the total energy, and gradients for two states

    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    energy_e : float
        Gaussian total calculated energy in Hartree for the excited state
    grad_e : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 etc. Hartree/Bohr for the excited state
    energy_g : float
        Gaussian total calculated energy in Hartree for the ground state
    grad_g : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 etc. Hartree/Bohr for the ground state

    """
    with open(in_name) as data:
        lines = data.readlines()
    grad_g = []
    grad_e = []
    for line in lines:
        if line.strip():
            if line.strip()[0].isalpha():
                reading_e = False
                reading_g = False
            if "( 1)     EIGENVALUE" in line:
                energy_g = float(line.split()[3])
            if "( 2)     EIGENVALUE" in line:
                energy_e = float(line.split()[3])
            if reading_g:
                for num in line.split():
                    grad_g.append(float(num))
            if reading_e:
                for num in line.split():
                    grad_e.append(float(num))
            if "Gradient of iOther State" == line.strip():
                reading_g = True
            if "Gradient of iVec State." == line.strip():
                reading_e = True

    grad_e = np.array(grad_e)
    grad_g = np.array(grad_g)
    return energy_e, grad_e, energy_g, grad_g


def read_g_dens(in_file, total_ci=False):
    """
    Read the density matrix from a gaussian fchk file

    Parameters
    ----------
    in_file : str
        Name of the file to read
    total_ci : bool optional
        If true, reads the total CI (excited state for e.g. TD-DFT) density.
        Otherwise it will only read the ground state density.
    Returns
    -------
    dens_mat : numpy matrix
        The density matrix

    """
    if total_ci:
        keyword = "Total CI Density"
    else:
        keyword = "Total SCF Density"

    with open(in_file) as to_read:
        reading = False
        entries = []
        for line in to_read:
            if line[0].isalpha():
                reading = False
            if reading == True:
                for num in [float(i) for i in line.split()]:
                    entries.append(num)
            if keyword in line:
                reading = True

    # there are (N^2+N)/2 entries in the half of an N x N matrix, the number of
    # the number of orbitals is therefore:
    n_orb = int((-1 + np.sqrt(1 + 8 * len(entries))) / 2)

    dens_mat = np.zeros((n_orb, n_orb))

    count = 0
    for i in range(n_orb):
        for j in range(i + 1):
            dens_mat[i][j] = entries[count]
            count += 1

    return dens_mat


def read_vectors(in_file):
    """
    Read a set of unit cell vectors from a formatted vector file

    Parameters
    ----------
    in_file : str
        Input file name
    Returns
    -------
    vectors : 3 x 3 numpy array
        Unit cell vectors where vector a is vectors[0], b is vectors[1], c is
        vectors[2]

    """
    vectors = np.loadtxt(in_file)
    if len(vectors) != 3:
        raise ValueError("The lattice vector file does not have 3 vectors")
    return vectors


def read_cube(in_file):
    """
    Read a cube file and return a Mol and a CubeGrid object

    Parameters
    ----------
    in_file : str
        Input file name
    Returns
    -------
    out_mol : Mol object
        The atoms in the cube file
    out_cub : CubeGrid object
        The grid in the cube file where all distances are in Angstrom

    """
    vectors = np.zeros((3, 3))
    xyz_nums = [0, 0, 0]
    values = []

    out_mol = Mol([])
    ind = 0
    natoms = 0
    with open(in_file) as lines:
        for line in lines:
            if ind == 2:
                natoms = int(line.split()[0])
                origin = np.array([float(i)
                                   for i in line.split()[1:]]) / pt.bohrconv
            if ind == 3:
                xyz_nums[0] = int(line.split()[0])
                vectors[0] = np.array([float(i)
                                       for i in line.split()[1:]]) / pt.bohrconv
            if ind == 4:
                xyz_nums[1] = int(line.split()[0])
                vectors[1] = np.array([float(i)
                                       for i in line.split()[1:]]) / pt.bohrconv
            if ind == 5:
                xyz_nums[2] = int(line.split()[0])
                vectors[2] = np.array([float(i)
                                       for i in line.split()[1:]]) / pt.bohrconv
                out_cub = CubeGrid(vectors, xyz_nums[0], xyz_nums[
                                   1], xyz_nums[2], origin)
                out_cub.set_grid_coord()
            if 6 <= ind < (6 + natoms):
                line_s = line.split()
                new_atom = Atom()
                new_atom.elem = per.num_to_elem(int(line_s[0]))
                new_atom.set_pos([float(i) / pt.bohrconv for i in line_s[2:]])
                out_mol.append(new_atom)
            if ind >= (6 + natoms):
                values.extend([float(i) for i in line.split()])
            ind += 1
    values_arr = np.array(values)
    out_cub.grid[:, 3] = values_arr
    return out_cub, out_mol

def read_dynamics(output_file):
    """
    This function is very similar to read_molcas. The difference is that in this functio
    all the states requested in the Molcas calculation are considered, whereas in read_molcas
    only the state selected is considered. I think that only one function can be written for
    both purposes (dynamics and optimization)
    """
    # Initialize variables
    energies = []
    gradients = np.array([])
    
    reading = False
    in_table = False

    with open(output_file,'r') as rf:
        rf_lines = rf.readlines()
    
    for line in rf_lines:
        if line.strip():
            # Energies
            if """::    RASSCF root number""" in line:
                energies.append(float(line.split()[-1]))
            # Gradients
            if reading:
                if len(line.split()) == 4:
                    if len(line.split()[0]) > 1:
                        if line.split()[0][1].isdigit() or line.split()[0][2].isdigit():
                            in_table = True
                            nums = [float(i) for i in line.split()[1:4]]
                            gradients = np.concatenate((gradients, nums))
            if in_table and "----------" in line:
                reading = False
                in_table = False
            if "Molecular gradients" in line:
#                print("Found molecular gradient output")
                reading = True
    energies = np.reshape(np.array(energies),(-1,1))

    return energies, gradients, energies[0]

def read_velocities(vel_file):
    """
    Read atomic velocities from formatted input file in [units]
    
    Parameters
    ----------
    vel_file : str
        Velocities file name
    Returns
    -------
    velocities : 3 x N Numpy array (N = number of atoms)
        Cartesian components of atomic velocities
    """
    velocities = []
    
    with open(vel_file,'r') as vf:
        vf_lines = vf.readlines()
        
        for line in vf_lines:
            (X,Y,Z) = line.split()
            velocities.append([float(X), float(Y), float(Z)])
    
    return velocities

def read_dyn_restart(restart_file):
    """
    Read gradients from formatted input file in atomic units
    The gradiends comes from the output of a trajectory that
    ended with error ("info_dyn_restart" file)
 
    Parameters
    ----------
    grad_file : str
        Gradients file name
    Returns
    -------
    Gradients : nstates x 3 N Numpy array (N = number of atoms)
        Cartesian components of atomic gradients

    """
    gradients = np.array([])

    reading = False
    in_table = False
    with open(restart_file,'r') as gf:
        gf_lines = gf.readlines()

    for line in gf_lines:
        if line.strip():
            if "Checkpoint at step" in line:
                curr_step = int(line.split()[-1])
            if "Total Energy at step 1:" in line:
                Eini = float(line.split()[-1])
            if "Gradients at previous step" in line:
                reading = True
        if reading:
            if len(line.split()) == 3:
                if len(line.split()[0]) > 1:
                    if line.split()[0][1].isdigit() or line.split()[0][2].isdigit:
                        in_table = True
                        nums = [float(i) for i in line.split()[:]]
                        gradients = np.concatenate((gradients,nums))
        if in_table and "----------" in line:
            reading= False
            in_table = False

    return curr_step, Eini, gradients  
         

def read_molcas_nacs(in_file):
    """
    Read Non adiabatic couplings from a Molcas output file with an arbitrary number of roots.
    The way that nacs are read for a case of N states is: St1 St2, St1 St3,.., 1St StN, 
    St2 St3,.., St2 StN, etc.
 
    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    nacs : numpy array of floats
        nacs in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    """
    # Initialize variables
    nacs = np.array([])

    reading = False
    in_table = False

    with open(in_file,'r') as rf:
        rf_lines = rf.readlines()

    for line in rf_lines:
        if line.strip():
            if "CI derivative coupling" in line:
#            if "CSF derivative coupling" in line:
                reading = True
            if reading:
                if len(line.split()) == 4:
                    if len(line.split()[0]) > 1:
                        if line.split()[0][1].isdigit() or line.split()[0][2].isdigit():
                            in_table = True
                            nums = [float(i) for i in line.split()[1:4]]
                            nacs = np.concatenate((nacs, nums))
            if in_table and "----------" in line:
                reading = False
                in_table = False

    return nacs

def read_gaussian_nacs(in_file):
    """
    Read Non adiabatic couplings between S0 and Sn from a Gaussian.fchk file
 
    Parameters
    ----------
    in_name : str
        Name of the file to read
    Returns
    -------
    nacs : numpy array of floats
        nacs in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    """
    # Initialize variables
    nacs = []
    reading = False

    with open(in_file,'r') as rf:
        rf_lines = rf.readlines()

    for line in rf_lines:
        if line[0].isalpha():
            reading = False
        if reading == True:
            for num in map(float, line.split()):
                nacs.append(num)
        if line.startswith("Nonadiabatic coupling"):
            reading = True
    nacs = np.array(nacs)
#    
    return nacs
#
def SF_dft_qchem(in_name):
    """
    Reads a qchem .out file.

    Returns the total energy, gradients and ground state energy
    for a spin-flip DFT calculation. 

    Parameters
    ----------
    in_name : str
        Name of the file to read

    Returns
    -------
    energy : float
        QChem total CDFT energy in Hartree
    grad : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 in Hartree/Bohr
    scf_energy : float
        QChem CDFT calculated energy in Hartree

    """
    with open(in_name) as qchem_file:
        content = qchem_file.readlines()

    grad = []
    x_grad = []
    y_grad = []
    z_grad = []

    reading=False
    reading_SASF = False
    read_grad = False
    count = 0

    for i, line in enumerate(content):

        # Reads the ground state SCF energy
        if line.strip().startswith("Total energy for") and line.split()[4] == "1:":
#        if line.split()[0] == "Total" and line.split()[4] == "1:":
        
#        if "Total energy for state  1:" in line:
            scf_energy = float(line.split()[5])
            continue

        # Reads the total energy from a SF-DFT output
        if "State Energy is" in line.strip() or "SASF ENERGY for OPT =" in line.strip():
#        if "State Energy is" in line.strip():
            total_energy = float(line.split()[-1])
            continue
        
        # Gradient starts after this line in DFT/CDFT
        # Set reading to True to start reading after this line
        if line.strip() == "Gradient of SCF Energy":
            reading = True
            continue

        # Gradient starts after this line in SF-DFT
        # Set reading to True after this line is found
        elif line.strip().startswith("Gradient of the state energy"):
            reading = True
            continue
     
        # Gradient starts after this line in SA-SF-DFT
        # Set reading to True after this line is found
        elif line.strip().startswith("SA-SF-DFT Gradient"):
            reading_SASF = True
            read_grad = True
            continue

        elif line.strip().startswith("SA-SF-DFT time:"):
            read_grad = False
        
        # Gradient output ends on the line before this line is found
        # Reading false (for TDDFT)
        elif line.strip().startswith("Gradient time:"):
            reading = False
            continue

        if reading_SASF == True and read_grad == True:
            if len(line.split()) == 4 and line.split()[0][-1].isdigit():
                #if len(line.split()[0]) > 1:
#                if line.split()[0][0].isdigit():
#                read_grad = True
#                    if line.split()[0][-1].isdigit():
#                if read_grad==True:
                print("FOUND read_grad TRUE")
                nums = [float(i) for i in line.split()[1:]]
                grad = np.concatenate((grad, nums))

        if reading:
            if count % 4 == 0:
                count += 1
                continue

            elif count == 1:
                for val in line.split()[1:]:
                    x_grad.append(float(val))
                count += 1
                continue

            elif count == 2:
                for val in line.split()[1:]:
                    y_grad.append(float(val))
                count += 1
                continue

            elif count == 3:
                for val in line.split()[1:]:
                    z_grad.append(float(val))
                count = 0
                continue

    if reading_SASF == False:
        for (x,y,z) in zip(x_grad, y_grad, z_grad):
            grad.append(x)
            grad.append(y)
            grad.append(z)

    grad = np.array(grad)
    
    print("total_energy=",total_energy)
    print("scf_energy=",scf_energy) 
    print("")
    print("GRAD")
    print(grad)

    return total_energy, grad, scf_energy

def read_qchem_out(in_name):
    """
    Reads a qchem .out file.

    Returns the total energy, gradients and ground state energy.

    Parameters
    ----------
    in_name : str
        Name of the file to read

    Returns
    -------
    energy : float
        QChem total CDFT energy in Hartree
    grad : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 in Hartree/Bohr
    scf_energy : float
        QChem CDFT calculated energy in Hartree

    """
    with open(in_name) as qchem_file:
        content = qchem_file.readlines()
        
    grad = []
    x_grad = []
    y_grad = []
    z_grad = []
    
    reading=False
    spin_flip = False
    count = 0
    true = ["true", "yes", "1"]
    
    
    for i, line in enumerate(content):

        #Check if it is a SF-DFT calculation and go to 
        #a different reading function in that case

        if line.strip().lower().startswith(("spin_flip","sasf_rpa")):
#        if line.strip().startswith("SPIN_FLIP") or line.strip().startswith("spin_flip"):
            spin_flip = str(line.split()[1])
            if spin_flip.lower() in true:
                spin_flip = True
                total_energy, grad, scf_energy = SF_dft_qchem(in_name)
                break

        # Reads the ground state SCF energy
        if line.strip().startswith("SCF   energy in the final basis set"):
            scf_energy = float(line.split()[8])
            continue
        
        # For the cases of ground state DFT and constrained DFT
        # This is necessary as fromage reads in the total energy to print
        # a gap, which in excited state calculations is the excitation 
        # energy for the chosen state
        # For CDFT/DFT, total energy is equal to SCF energy and the
        # resulting gap is 0
        
        if line.strip().startswith("Total energy in the final basis set"):
            total_energy = float(line.split()[8])
            continue
        
        # Reads the total energy from a TDDFT output
        # Where RPA is True
        if "RPA" and "State Energy is" in line.strip():
            total_energy = float(line.split()[5])
            continue
        
        # Gradient starts after this line in DFT/CDFT
        # Set reading to True to start reading after this line
        if line.strip() == "Gradient of SCF Energy":
            reading = True
            continue
            
        # Gradient starts after this line in TDDFT
        # Set reading to True after this line is found
        elif line.strip().startswith("Gradient of the state energy"):
            reading = True
            continue
        
        # Gradient output ends on the line before this line is found
        # Reading false (for DFT/CDFT)
        if line.strip().startswith("Max gradient component"):
            # Turn off reading as the gradient should be fully 
            # read in by this line
            reading = False
            continue
        
        # Gradient output ends on the line before this line is found
        # Reading false (for TDDFT)
        elif line.strip().startswith("Gradient time:"):
            reading = False
            continue
        
        # The following blocks of code are the same regardless of chosen
        # method
        # Qchem also has a consistent format of printing gradients across
        # its methods
        # Units are hartree/bohr (atomic units)
        # First line is atom number (order in input file)
        # 2nd-4th line are x y and z components of gradient
        #      1 2 3 4 5 6
        #      x x x x x x
        #      y y y y y y 
        #      z z z z z z 
        # Largest number of atoms in a single row is 6 
        # Then the gradients for the next row of atoms are printed
        
        if reading:
            if count % 4 == 0:
                count += 1 
                continue
            
            elif count == 1:
                for val in line.split()[1:]:
                    x_grad.append(float(val))
                count += 1
                continue
                
            elif count == 2:
                for val in line.split()[1:]:
                    y_grad.append(float(val))
                count += 1
                continue
                
            elif count == 3:
                for val in line.split()[1:]:
                    z_grad.append(float(val))
                count = 0
                continue
                      
    if spin_flip == False:
        for (x,y,z) in zip(x_grad, y_grad, z_grad):
            grad.append(x)
            grad.append(y)
            grad.append(z)
        
        grad = np.array(grad) 
    
    return total_energy, grad, scf_energy


def read_nwchem_DFT_out(in_name):
    """
    Reads DFT, CDFT and TD-DFT output files from NWChem.

    N.B Gradients are not implemented for CDFT spin constraints, only charge constraints.

    Returns the total energy, gradients and ground state energy.

    Parameters
    ----------
    in_name : str
        Name of the file to read

    Returns
    -------
    energy : float
        NWChem total energy in Hartree
    grad : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 in Hartree/Bohr
    scf_energy : float
        NWChem ground state SCF energy in Hartree
    """
    with open(in_name) as nwchem_file:
        content = nwchem_file.readlines()

    grad = []

    reading = False
    total_energy = None

    for i, line in enumerate(content):
        # Reads the ground state SCF energy or CDFT energy
        if line.strip().startswith("Total DFT"):
            scf_energy = float(line.split()[4])

        # Reads the excited state energy for the chosen state
        if line.strip().startswith("Excited state energy"):
            total_energy = float(line.split()[4])

        # Gradient starts after this line in DFT/CDFT/TDDFT
        # Set reading to True to start reading after this line
        if "ENERGY GRADIENTS" in line.strip():
            reading = True

        # Gradient format is as follows (in atomic units)
        # atomnumber elem coordinates gradients
        # e.g.
        # 1 C x1_coord y1_coord z1_coord x1_grad y1_grad z1_grad
        # 2 O x2_coord y2_coord z2_coord x2_grad y2_grad z2_grad

        # check that line corresponds to gradient
        if reading:
            # grad line has 8 elements
            if len(line.split()) == 8:
                try:
                    grad_vals = [float(i) for i in line.split()[5:]]
                    grad = grad + grad_vals
                # ValueError means we are past the gradient output
                # so set reading to False
                except ValueError:
                    reading = False
                    break

    # For CDFT/DFT we need to set total energy equal to SCF energy
    # This is necessary as fromage reads in the total energy to print
    # a gap, which in excited state calculations is the excitation
    # energy for the chosen state. In CDFT/DFT the gap is 0.

    if total_energy == None:
        total_energy = scf_energy

    grad = np.array(grad)

    return total_energy, grad, scf_energy

def mopac_fomo_ci_out(in_name):
    """
    Reads FOMO-CI output files from MOPAC

    Returns the total energy, gradients and ground state energy.

    Parameters
    ----------
    in_name : str
        Name of the file to read

    Returns
    -------
    state_energy : float
        MOPAC-FOMO-CI total energy in Hartree
    grad : list of floats
        The gradients in form x1,y1,z1,x2,y2,z2 in Hartree/Bohr
    gr_energy : float
        MOPAC-FOMO-CI ground state SCF energy in Hartree
    """
    with open(in_name) as fomo_ci_file:
        content = fomo_ci_file.readlines()

    grad = []
    reading = False
    state_energy = None

    for i, line in enumerate(content):

        # Reads the excited state energy for the chosen state in a.u.
        if "KINETIC ENERGY=" in line:
            state_energy = float(line.split()[-2])
        # Gradients of the selected state in a.u. (Hartree / bohr)
        if "GRADIENT OF THE POTENTIAL" in line:
            reading = True
            continue
        if reading:
            if len(line.split()) == 3:
                nums = [float(i) for i in line.split()[:]]
                grad = np.concatenate((grad,nums))
            else:
                reading = False
        # GS scf energy. NB, if the multiplicity is set to TRIPLETS,
        # the energy here corresponds to T1, not to S0.
        if "CI VECTORS" in line:
            gr_energy = float(content[i+5].split()[0])

    return state_energy, grad, gr_energy
#
def read_orca_out(in_name):
    """
    Read energies and gradients from an Orca output

    Parameters
    ----------
    in_name : str
        name of the file to read
    Returns
    -------
    state_energy : float
        The energy (in Hartrees) of the required sate
    grad : numpy array of floats
        Energy gradients in the form x1,y1,z1,x2,y2,z2 etc. in Hartree/Bohr
    gr_energy : float
        Ground state energy in Hartree

    """
    with open(in_name) as Orca_file:
        content = Orca_file.readlines()

    grad = []
    gr_energy=None
    state_energy=None

    for i, line in enumerate(content):
        if "Total Energy" in line:
            gr_energy = float(line.split()[3])
        if "FINAL SINGLE POINT ENERGY" in line:
            state_energy = float(line.split()[4])
            # For a ground state calculation:
            # ex_energy = gr_energy
        if "CARTESIAN GRADIENT" in line:
            orig_line = i
            break
    for line in content[orig_line + 3:]:
        if len(line.split()) == 6:
            atom_grads = [float(i) for i in line.split()[3:]]
            grad = np.concatenate((grad,atom_grads))
        else:
            break
#
    if gr_energy == None:
        gr_energy = state_energy
#
    return state_energy, grad, gr_energy
###################################################################
#####################  Read Hessians  #############################
"""
 The only input is the name of the file where the Hessian is stored

"""
def read_hessian_xtb(in_name):
    with open(in_name) as data:
        lines = data.readlines()
    hess = []
    for line in lines[1:]:
        for num in map(float, line.split()):
            hess.append(num)
    hess_dim = int(np.sqrt(len(hess)))
    hess = np.array(hess).reshape(hess_dim,hess_dim)
    return hess

def read_hessian_dftb(in_name):
    """
    The hessian in dftb+ is presented as:
    d^2E/(dx1dx1) d^2E/(dy1dx1) d^2E/(dz1dx1) 
    d^2E/(dx2dx1) d^2E/(dy2dx1) d^2E/(dz2dx1)
    .
    .
    d^2E/(dxNdx1) d^2E/(dyNdx1) d^2E/(dzNdx1)
    .
    .
    d^2E/(dxNdxN) d^2E/(dyNdxN) d^2E/(dzNdxN)

    """
    with open(in_name) as data:
        lines = data.readlines()
    hess = []
    for line in lines:
        for num in map(float, line.split()):
            hess.append(num)
    hess_dim = int(np.sqrt(len(hess)))
    hess = np.array(hess).reshape(hess_dim,hess_dim)
    return hess

def read_hessian_g_fchk(in_name):
    with open(in_name) as data:
        lines = data.readlines()
    hess_lt = []
    reading = False
    for line in lines:
        if line[0].isalpha():
            reading = False
        if reading == True:
            # Get the lower triangle of the Hessian matrix
            for num in map(float, line.split()):
                hess_lt.append(num)
        if line.startswith("Cartesian Force Constants"):
            reading = True
        if line.startswith("Number of atoms"):
            natoms = int(line.split()[-1])
            hess_dim = int(3.* natoms)
    hess_lt = np.array(hess_lt)
    hess = np.zeros((hess_dim,hess_dim))
    # convert the lower triangle matrix to the symetric full Hessian
    hess[np.tril_indices(hess.shape[0], k = 0)] = hess_lt
    return hess
#
def read_hessian_turbomole(in_name):
    with open(in_name) as data:
        lines = data.readlines()
    hess = []
    for line in lines[1:-1]:
        for num in map(float, line.split()[2:]):
            hess.append(num)
    hess_dim = int(np.sqrt(len(hess)))
    hess = np.array(hess).reshape(hess_dim,hess_dim)
    return hess
#
def read_hessian_molcas(in_name):
    """
    In Molcas, the hessian is obtained from a binary h5 file. Hence,
    the module h5py must be imported
    """
    import h5py
    f = h5py.File(in_name, 'r')
    nuclei_key = list(f.keys())[0]
    hessian_key = list(f.keys())[-1]
    natoms = len(list(f[nuclei_key]))
    # Get the lower triangle of the Hessian matrix
    hess_lt = list(f[hessian_key])
    f.close()
    hess_dim = int(3.* natoms)
    hess_lt = np.array(hess_lt)
    hess = np.zeros((hess_dim,hess_dim))
    # convert the lower triangle matrix to the symetric full Hessian
    hess[np.tril_indices(hess.shape[0], k = 0)] = hess_lt
    return hess
#
def read_hessian_orca(in_name):
    """
    Read the hessian matrix computed by Orca
    """
    with open(in_name, 'r') as file:
        lines = file.readlines()

    # Find the line index where $hessian is located
#    hessian_index = -1
#    for i, line in enumerate(lines):
#        if "$hessian" in line:
#            hessian_index = i
#            break

#    if hessian_index == -1:
#        raise ValueError(f"Hessian data not found in the file {file_path}")

    hessian_index = next((i for i, line in enumerate(lines) if "$hessian" in line), None)
    if hessian_index is None:
        raise ValueError(f"Hessian data not found in the file {in_name}")

    dim = int(lines[hessian_index + 1].strip())
    hess = np.zeros((dim, dim))

    current_row = 0
    total_col_index = 0  # Total column index across all blocks

    for line in lines[hessian_index + 2:]:
        if line.strip(): 
            values = [float(val) for val in line.split()[1:]]

            for i, value in enumerate(values):
                if total_col_index + i < dim:
                    hess[current_row, total_col_index + i] = value

            current_row += 1

            if current_row == dim:
                total_col_index += len(values)
                current_row = 0

                if total_col_index >= dim:
                    break

    return hess


"""
 To implement: Q-Chem - NWChem - dftbplus
"""
###################################################################
################  Read dipole derivatives  ########################
def read_dfrb_mu(in_name):
    """
    Read the dipole derivatives computed with dftb+
    """
    with open(in_name, 'r') as data:
        lines = data.readlines()
    d_mu = []
    dim = len(lines)
    for line in lines:
        for num in map(float, line.split()):
            d_mu.append(num)
    d_mu = np.array(hess).reshape(dim,3)

    return d_mu

def read_gauss_mu(in_name):
    """
    Read the dipole derivatives computed with G16/G09
    """
    start_reading = False
    nums = []
    dim = 0

    for line in in_name.splitlines():
        if 'Dipole Derivatives' in line:
            start_reading = True
            dim = int(line.split()[-1])
            continue
        if 'Polarizability' in line:
            break
        if start_reading:
            nums.extend([float(num) for num in line.split()])

    return np.array(numbers).reshape((int(dim / 3), 3))

def read_turbo_mu(in_name):
    """
    Read the dipole derivatives computed with Turbomole
    """

    with open(in_name) as data:
        lines = data.readlines()

    d_mu = []
    for line in lines[1:-1]:
        for num in map(float, line.split()[2:]):
            d_mu.append(num)
    dim = int(len(lines) - 2)

    return np.array(d_mu).reshape(dim,3)

def read_orca_mu(in_name):
    """
    Read the dipole derivatives computed with Turbomole
    """

    lines = in_name.splitlines()
    numbers = []
    dim = 0
    read_data = False

    for i, line in enumerate(lines):
        if '$dipole_derivatives' in line:
            dim = int(lines[i + 1].strip())
            read_data = True
            continue
        if read_data:
            nums.extend([float(num) for num in line.split()])
            if len(nums) // 3 == dim:
                break

    return np.array(numbers).reshape((dim, 3))
