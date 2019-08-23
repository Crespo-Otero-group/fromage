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
    trajectory.

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

    return energy, grad


def read_molcas(in_name):
    """
    Read energies and gradients from a Molcas .log file with 2 roots.

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
    with open(in_name) as data:
        lines = data.readlines()

    grad = np.array([])
    ex_energy = None
    gr_energy = None

    reading = False

    for line in lines:
        if line.strip():
            # Energies
            if "RASSCF root number  1 Total energy:" in line:
                gr_energy = float(line.split()[-1])
            if "RASSCF root number  2 Total energy:" in line:
                ex_energy = float(line.split()[-1])
            # Gradients
            if "Molecular gradients" in line:
                reading = True
            if reading:
                if len(line.split()) == 4 and line.split()[0][0].isalpha():
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
    ex_energy : None
        There are no excited state gradients in dftb+ so this will always be empty
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
                atom_grads = [float(i) for i in line.split()]
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
