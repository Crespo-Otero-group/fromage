"""Functions for creating files needed for other software"""

import numpy as np
import subprocess
import os
from random import randint

def write_cp2k(in_name, file_name, vectors, atoms, temp_name):
    """
    Make a cp2k input file from a template.

    To use this, first make a template file with the name as XXX__NAME__XXX,
    the vectors as XXX__AVEC__XXX, XXX__BVEC__XXX and XXX__CVEC__XXX and the
    positions as XXX__POS__XXX.

    Parameters
    ----------
    in_name : string
        Project name
    file_name : string
        Name of the cp2k input file
    vectors : 3x3 matrix
        Lattice vectors
    atoms : list of Atom objects
        Atoms in the system
    temp_name : string
        Name of the template file


    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    # strings for each lattice vector
    aVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[0])
    bVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[1])
    cVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[2])

    cp2k_in = open(file_name, "w")

    for line in temp_content:

        # writes the name of the calculation at the top of the file
        if "XXX__NAME__XXX" in line:
            cp2k_in.write(line.replace("XXX__NAME__XXX", in_name))

        # replace the tags with the coordinates of lattice vectors
        elif "XXX__AVEC__XXX" in line:
            cp2k_in.write(line.replace("XXX__AVEC__XXX", aVec))
        elif "XXX__BVEC__XXX" in line:
            cp2k_in.write(line.replace("XXX__BVEC__XXX", bVec))
        elif "XXX__CVEC__XXX" in line:
            cp2k_in.write(line.replace("XXX__CVEC__XXX", cVec))

        # writes atomic coordinates
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                line_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z)
                cp2k_in.write(line_str + "\n")

        else:  # if no tag is found
            cp2k_in.write(line)

    cp2k_in.close()
    return


def write_xyz(in_name, atoms, char=False):
    """
    Write an xyz file.

    Parameters
    ----------
    in_name : string
        Name of the xyz file. Include the file extension, e.g. "molecule.xyz"
    atoms : list of Atom objects
        Atms to write
    char : optional bool
        Write the charge of the atom in the 4th column

    """
    out_file = open(in_name, "w")
    out_file.write(str(len(atoms)) + "\n")
    out_file.write(in_name + "\n")

    for atom in atoms:
        if char:
            out_file.write(str(atom) + "\n")
        else:
            out_file.write(atom.xyz_str() + "\n")
    out_file.close()
    return


def write_uc(in_name, vectors, aN, bN, cN, atoms):
    """
    Write a .uc file for Ewald.c.

    .uc files are written in fractional coordinates instead of Cartesian

    Parameters
    ----------
    in_name : string
        Name of the file to be written
    vectors : 3x3 matrix
        Lattice vectors
    aN,bN,cN : ints
        Number of times each lattice vector should be multiplied
    atoms : list of Atom objects
        Unit cell atoms for the file

    """
    line1 = vectors[0].tolist() + [aN]
    line2 = vectors[1].tolist() + [bN]
    line3 = vectors[2].tolist() + [cN]
    out_file = open(in_name, "w")
    out_file.write("{:10.6f} {:10.6f} {:10.6f} {:10d}".format(*line1) + "\n")
    out_file.write("{:10.6f} {:10.6f} {:10.6f} {:10d}".format(*line2) + "\n")
    out_file.write("{:10.6f} {:10.6f} {:10.6f} {:10d}".format(*line3) + "\n")

    # transpose to get the transformation matrix
    M = np.transpose(vectors)
    # inverse transformation matrix
    U = np.linalg.inv(M)

    for atom in atoms:
        # change of basis transformation
        dir_pos = [atom.x, atom.y, atom.z]
        frac_pos = np.dot(U, dir_pos).tolist()
        for coord in frac_pos:
            # if the coordinate is out of range
            if coord < 0 or coord > 1:
                # translate it to the range [0,1]
                frac_pos[frac_pos.index(coord)] = coord % 1
        str_line = "{:10.6f} {:10.6f} {:10.6f} {:14.10f} {:>6}".format(
            *frac_pos + [atom.q] + [atom.elem]) + "\n"
        out_file.write(str_line)
    out_file.close()
    return


def write_qc(in_name, atoms):
    "Write a .qc file for Ewald.c."
    out_file = open(in_name, "w")
    for atom in atoms:
        str_line = "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
            atom.elem, atom.x, atom.y, atom.z, atom.q) + "\n"
        out_file.write(str_line)
    out_file.close()
    return


def write_ew_in(in_name, file_name, nChk, nAt):
    """
    Write an input file for Ewald.c.

    Remember to include and pre or post fixes to the file name. Following the
    convention from the paper, it would be something like "in.ewald.somename"

    Parameters
    ----------
    in_name : string
        Name of the project
    file_name : string
        Name of the input file to write
    nChk : int
        Number of random checkpoints in zone 1
    nAt : int
        Number of fixed charge atoms in zone 1+2

    """
    out_file = open((file_name), "w")
    out_file.write(in_name + "\n")
    out_file.write(str(nChk) + "\n")
    out_file.write(str(nAt) + "\n")
    out_file.write("0\n")
    out_file.close()
    return


def write_seed():
    """Write a seedfile for Ewald.c"""
    out_file = open("seedfile", "w")
    seed1 = randint(1, 2**31 - 86)
    seed2 = randint(1, 2**31 - 250)
    out_file.write(str(seed1) + " " + str(seed2))
    out_file.close()


def write_gauss(file_name, atoms, points, temp_name, proj_name='gaussian', freq=None, state=None, states=None):
    """
    Write a Gaussian input file.

    A template file needs to be prepared which has the name as XXX__NAME__XXX,
    the positions as XXX__POS__XXX and the (optional) charges as
    XXX__CHARGES__XXX.

    Parameters
    ----------
    file_name : str
        Name of the Gaussian input file to be written
    atoms : list of Atom objects
        Atoms to be calculated with Gaussian
    points : list of Atom objects or None
        The Atom objects should have a charge. If there is no XXX__CHARGES__XXX
        in the template file, this doesn't matter and can be None
    temp_name : str
        Name of the template file
    proj_name : str
        Project name, default 'gaussian'

    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")

    for line in temp_content:
        if "XXX__NAME__XXX" in line:
            out_file.write(line.replace("XXX__NAME__XXX", proj_name))
        if line.startswith("#"):
            if freq is not None:
                if "force" in line.strip():
                    out_file.write(line.replace("force", "freq"))
                elif "Force" in line.strip():
                    out_file.write(line.replace("Force", "freq"))                
                else:
                    line = line.strip() + ' freq\n'
                    out_file.write(line)
            if not ("symmetry=none" in line or "Nosymm" in line or "nosymm" in line):
                line += " symmetry=none"
            if "&STATE" in line and state != None:
                if "&NSTATES" in line and states != None:
                    curr_state = '%s' % (state - 1)
                    nstates = '%s' % (int(np.sum(states))-1)
                    modified_line = line.replace("&STATE", curr_state).replace("&NSTATES", nstates)
                    out_file.write(modified_line)
            else:
                out_file.write(line)
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        elif "XXX__CHARGES__XXX" in line:
            for point in points:
                point_str = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    point.x, point.y, point.z, point.q) + "\n"
                out_file.write(point_str)
        else:
            out_file.write(line)
    out_file.close()
    return

def write_dftb(file_name, atoms, points, temp_name, proj_name='dftb',freq=None):
    """
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")

    for line in temp_content:
        if "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = " {:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        else:
            out_file.write(line)
    out_file.close()

    if freq is not None:
#        subprocess.call("mv dftb_in.hsd dftb_in.hsd.temp", shell=True)
        write_dftb_freq_calc('dftb_in.hsd')

    return

def write_dftb_freq_calc(file_name):
    """
    Updates the 'Driver' section in the dftb+ input file.

    Parameters:
    file_path (str): Path to the file to be modified.
    """
    with open(file_name, 'r') as file:
        lines = file.readlines()

    if any('Driver = SecondDerivatives{' in line for line in lines):
        return

    new_driver_block = [
        "Driver = SecondDerivatives {\n",
        "    Atoms = 1:-1\n",
        "    Delta = 1e-5\n",
        "    }\n"
    ]

    with open(file_name, 'w') as file:
        for line in lines:
            if line.strip() == 'Driver = {}':
                file.writelines(new_driver_block)
            else:
                file.write(line)

    return

def write_dftb_charges(file_name, points):
    """
    """
    out_file = open(file_name, "w")

    for point in points:
        point_str = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
            point.x, point.y, point.z, point.q) + "\n"
        out_file.write(point_str)
    out_file.close()

    return

def write_xtb(file_name, atoms, points, temp_name, proj_name='xtb'):
    """
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")

    for line in temp_content:
        if "XXX__NAME__XXX" in line:
            out_file.write(line.replace("XXX__NAME__XXX", proj_name))
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        else:
            out_file.write(line)
    out_file.close()
    if len(points) > 0:
        out_file = open("xtb_charge.pc","w")
        out_file.write("%s \n" % len(points))
        for point in points:
            point_str = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                point.q, point.x, point.y, point.z) + "\n"
            out_file.write(point_str)
        out_file.close()
         
    return

def write_turbo_dyn(file_name,
                    temp_name,
                    state : int,
                    states : list,
                    singlestate : int,
                    nac_coupling : list,
                    soc_coupling : list):

    """
     Modify the Turbomole control file to include the calculation of the gradients for all
     the states considered in the dynamics or only for the current state, according to the
     user's preference

     All the input files for the force calculation with Turbomole are required. The  

    """

    with open(temp_name, 'r') as temp_file:
#        temp_content = temp_file.read().split('&')
        temp_content = temp_file.readlines()

    with open(file_name, 'w') as out_file:
        for line in temp_content:
            stripped_line = line.strip()  # Remove leading and trailing whitespace
  
            # Block for ADC2/CC2 calculations
            if '&NEXC' in stripped_line:
                nstates = '%s' % (int(np.sum(states))-1)
                nstart = '  nstart=%s' % (state - 1)
                npre = '  npre=%s' % (state - 1)
                nexc = '  nexc=%s' % (nstates)
                line = line.replace("&NEXC", nexc).replace("&NPRE", npre).replace("&NSTART", nstart)
            if '&GRAD' in stripped_line:
                if singlestate == 1:
                    gradient='  xgrad states(a %s)' % (state - 1)
                    line = line.replace('&GRAD', gradient)
                else:
                    gradient= 'xgrad states=all'
                    line = line.replace('&GRAD', gradient)

            # Block for TD-DFT calculations
            if '&SOES' in stripped_line:
                nstates = '%s' % (int(np.sum(states))-1)
                soes = '$soes all %s' % nstates
                line = line.replace('&SOES', soes)

            if '&EXOPT' in stripped_line:
                exopt = '$exopt %s' % (state - 1)
                line = line.replace('&EXOPT', exopt)

            out_file.write(line)

    return   


def write_molcas(file_name, temp_name, point_charges=None, freq=None):
    """
    This function is to update the point charges in the molcas
    input file when the QM' relaxation option is ON
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")

    for line in temp_content:
        if "XXX__NAME__XXX" in line:
            out_file.write(line.replace("XXX__NAME__XXX", proj_name))
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:10.6f} {:10.6f} {:10.6f} {:>6}".format(
                    atom.x, atom.y, atom.z, atom.elem) + "\n"
                out_file.write(atomStr)
        elif "XXX__CHARGES__XXX" in line:
            out_file.write("%s Angstrom \n" % len(point_charges))
            for point in point_charges:
                point_str = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    point.x, point.y, point.z, point.q) + "\n"
                out_file.write(point_str)
        else:
            out_file.write(line)
    out_file.close()
    
    if freq is not None:
        out_file = open(file_name, "a")
        out_file.write(" " + "\n")
        out_file.write("&MCKINLEY" + "\n")
        out_file.write("Perturbation = Hessian" + "\n")
        out_file.write("ShowHessian" + "\n")
        out_file.close()

    return

def write_molcas_free(file_name, 
                      temp_name,
                      state : int,
                      states : list,
                      singlestate : int,
                      nac_coupling : list,
                      soc_coupling : list,
                      point_flex : list):
    """
    Write a free format Molcas input file for dynamics run

    A tempalte file needs to be prepared which has the placeholder &GRAD and &NAC
    it may contain multiple sets of &RASSCF, &GRAD, and &NAC sections for different spin states

    -- by Jingbai Li 2022-02-22 and Federico Hernandez 2022-10-27
    
    Parameters
    ----------
    file_name : str
        Name of the Molcas input file to be written, default molcas.input
    temp_name : str
        Name of the template file ("molcas.temp")
    state : int
        Root number of the electronic state for which gradients should
        be followed for dynamics. Roots start with 1, so, for example,
        state = 2 is the first excited state.
    states : list
        List of number of states per spin multiplicity
    singlestate : int
        Flag to only compute gradient of the current state, the others will be zero to keep the (nstates, natom, 3) shape
    nac_coupling : list
        List of nac state pairs
    soc_coupling : list
        List of soc state pairs

    """

    with open(temp_name, 'r') as temp_file:
        temp_content = temp_file.read().split('&')

    # prepare grad section
    grad = []  # grad sections
    sect = []  # section index for each state
    indx = 0   # state index
    for s, ns in enumerate(states):
        sub = [] # subsections of grad
        for n in range(ns):
            indx += 1
            if singlestate == 1 and indx != state: # skip other state if only single state grad is requested
                alaska = ''
            else:
                alaska = 'ALASKA\nROOT=%s\n' % (n + 1)
#                alaska = 'ALASKA\nPNEW\nROOT=%s\n' % (n + 1)
            sub.append(alaska)
            sect.append(s)
        grad.append(sub)

    # prepare nac section
    nac = [[] for x in grad] # nac should have the same number of section as the grad
    if len(nac_coupling) > 0:
        for pair in nac_coupling:
            s1, s2 = pair   # two states
            alaska = 'ALASKA\nNAC=%s %s\n' % (s1 + 1, s2 + 1)
#            alaska = 'ALASKA\nPNEW\nNAC=%s %s\n' % (s1 + 1, s2 + 1)
            nac[sect[s1 - 1]].append(alaska)

    # prepare soc section
    soc = ['>>COPY   $WorkDir/$Project.JobIph   $WorkDir/JOB001\n', '>>COPY   $WorkDir/$Project.JobIph   $WorkDir/JOB002\n', '']
    if len(soc_coupling) > 0:
        na = states[0] # number of spin state a
        nb = states[1] # number of spin state b
        sa = [str(x + 1) for x in range(na)] # states of spin a
        sb = [str(x + 1) for x in range(nb)] # states of spin b
        sa = ' '.join(sa)
        sb = ' '.join(sb)
        soc[2] = 'RASSI\nNrofJobIph=2 %s %s;%s;%s\nSpinOrbit\nEJob\nSOCOupling=0\n' % (na, nb, sa, sb)

    # combine input
    input = []
    section = -1
    for n, line in enumerate(temp_content):
        if 'RASSCF' in line.upper():
            section += 1
        if 'GRAD' in line:
            for x in grad[section]:
                input.append(x)
        elif 'NAC' in line:
            for x in nac[section]:
                input.append(x)
        elif 'SOC' in line:
            input.append(soc[section])
        else:
            input.append(line)

    # add soc input, if soc_coupling != [], otherwise soc[2] == ''. 
    input.append(soc[2])

    # remove empty section and save
    input = '&'+'&'.join([x for x in input if x])
    input = input.replace('&>>', '>>') # remove an extra '&' before '>> COPY'

    with open(file_name, 'w') as out_file:
        out_file.write(input)

    return


def write_dynamics(file_name, temp_name, state : int, nstates : int):
    """
    Write a Molcas input file for a dynamics run.
    
    A tempalte file needs to be prepared which has the parameter Rlxroot
    set to XXX__STATE__XXX.
    
    Parameters
    ----------
    file_name : str
        Name of the Molcas input file to be written, default molcas.input
    temp_name : str
        Name of the template file
    state : int
        Root number of the electronic state for which gradients should
        be computed. Roots start with 1, so, for example, state = 2 is 
        the first excited state.
    
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, 'w')

    for line in temp_content:
        if "XXX__STATE__XXX" in line:
            out_file.write(line.replace("XXX__STATE__XXX", str(state)))
        elif "XXX__NSTATES__XXX" in line:
            out_file.write(line.replace("XXX__NSTATES__XXX", str(nstates)))
        else:
            out_file.write(line)
    out_file.close()
    return

def write_g_temp(file_name, fixed_atoms, points, temp_name, proj_name='gaussian'):
    """
    Write a Gaussian input template file.

    This serves to generate templates for write_gauss. You need a .template file
    which will generate a .temp file to be used in calculation. In this case a XXX__POS__XXX tag
    should be included in .template and won't be overwritten.

    Parameters
    ----------
    in_name : str
        Project name
    file_name : str
        Name of the Gaussian input file to be written
    fixed_atoms : list of Atom objects
        Atoms to be calculated with Gaussian with the -1 tag which makes them
        fixed in space
    points : list of Atom objects or None
        The Atom objects should have a charge. If there is no XXX__CHARGES__XXX
        in the template file, this doesn't matter and can be None
    temp_name : str
        Name of the template file
    pro_name : str
        Project name. Default 'gaussian'


    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")

    for line in temp_content:
        if "XXX__NAME__XXX" in line:
            out_file.write(line.replace("XXX__NAME__XXX", proj_name))
        elif "XXX__FIX__XXX" in line:
            for atom in fixed_atoms:
                atomStr = "{:>6} -1 {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        elif "XXX__CHARGES__XXX" in line:
            for point in points:
                point_str = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(
                    point.x, point.y, point.z, point.q) + "\n"
                out_file.write(point_str)
        else:
            out_file.write(line)
    out_file.close()
    return

def write_xtb_temp(file_name, fixed_atoms, points, temp_name, proj_name='xtb'):
    """
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, 'w')

    for line in temp_content:
        if "XXX__FIX__XXX" in line:
            for atom in fixed_atoms:
                atom_string = "{:10.6f} {:10.6f} {:10.6f} {:>6}\n".format(
                    atom.x, atom.y, atom.z, atom.elem)
                out_file.write(atom_string)
        else:
            out_file.write(line)
    out_file.close()

    charge_file = open("xtb_charge.pc", "w")

    charge_file.write(str(len(points)) + "\n")
    for point in points:
        point_str = "{:10.6f} {:10.6f} {:10.6f} {:10.6f}\n".format(
            point.q, point.x, point.y, point.z)
        charge_file.write(point_str)
    charge_file.close()

    return

def edit_vasp_pos(in_name, atoms):
    """
    Overwrite vasp POSCAR file.

    Parameters
    ----------
    in_name : str
        Name of the VASP file
    atoms : list of Atom objects
        Atoms to be written in the POSCAR file

    """
    with open(in_name + ".vasp") as vaspFile:
        content = vaspFile.readlines()

    out_file = open(in_name + ".new.vasp", "w")

    for line in content:
        out_file.write(line)
        if "Cartesian" in line:
            break
    for atom in atoms:
        atomStr = "{:10.6f} {:10.6f} {:10.6f}".format(
            atom.x, atom.y, atom.z) + "\n"
        out_file.write(atomStr)
    return


def write_qe(in_name, file_name, vectors, atoms, temp_name):
    """
    Write a Quantum Espresso (QE) input file.

    Parameters
    ----------
    in_name : str
        Project name
    file_name : str
        Name of the QE input file
    vectors : 3x3 matrix
        Lattice vectors
    atoms : list of Atom objects
        Atoms to be calculated with QE
    temp_name : str
        Name of the template file

    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    # strings for each lattice vector
    aVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[0])
    bVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[1])
    cVec = "{:10.6f} {:10.6f} {:10.6f}".format(*vectors[2])

    qe_in = open(file_name, "w")

    for line in temp_content:

        # writes the name of the calculation at the top of the file
        if "XXX__NAME__XXX" in line:
            qe_in.write(line.replace("XXX__NAME__XXX", in_name))

        # replace the tags with the coordinates of lattice vectors
        elif "XXX__AVEC__XXX" in line:
            qe_in.write(line.replace("XXX__AVEC__XXX", aVec))
        elif "XXX__BVEC__XXX" in line:
            qe_in.write(line.replace("XXX__BVEC__XXX", bVec))
        elif "XXX__CVEC__XXX" in line:
            qe_in.write(line.replace("XXX__CVEC__XXX", cVec))

        # writes atomic coordinates
        elif "XXX__POS__XXX" in line:
            for atom in atoms:
                line_str = "{:<6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z)
                qe_in.write(line_str + "\n")

        else:  # if no tag is found
            qe_in.write(line)
    qe_in.close()
    return


def write_pp(in_name, file_name, temp_name):
    """
    Write a Quantum Espresso .pp file

    This is a file for PP.x calculations which generate things like cube files

    Parameters
    ----------
    in_name : str
        Project name
    file_name : str
        Name of the .pp file. Include the extension e.g. "example.pp"
    temp_name :
        Name of the template file

    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    pp_in = open(file_name, "w")

    for line in temp_content:

        if "XXX__NAME__XXX" in line:
            pp_in.write(line.replace("XXX__NAME__XXX", in_name))
        else:
            pp_in.write(line)
    pp_in.close()
    return

def write_qchem(file_name, atoms, temp_name):
    """
    Write a QChem input file.

    A template file needs to be prepared which has XXX__POS__XXX
    instead of the atomic positions in the $molecule block.

    All other options also need to be specified in this file, e.g.
    symmetry false, and the force keyword to get the gradients.

    Parameters
    ----------
    file_name : str
        Name of the QChem input file to be written
    atoms : list of Atom objects
        Atoms to be calculated with QChem
    temp_name : str
        Name of the template file
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")
    # replaces XXX__POS__XXX with atom coordinates
    for line in temp_content:
        if "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        else:
            out_file.write(line)
    out_file.close()
    return

def write_nwchem(file_name, atoms, temp_name):
    """
    Write a NWChem input file. Currently, this function is redundant
    with the write_qchem function.

    A template file needs to be prepared which has XXX__POS__XXX
    instead of the atomic positions in the geometry.

    All other options must be specified in this file, e.g.
    symmetry, and the gradient keyword to get the gradients.

    Parameters
    ----------
    file_name : str
        Name of the NWChem input file to be written
    atoms : list of Atom objects
        Atoms to be calculated with NWChem
    temp_name : str
        Name of the template file
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")
    # replaces XXX__POS__XXX with atom coordinates
    for line in temp_content:
        if "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        else:
            out_file.write(line)
    out_file.close()
    return

def write_orca(file_name, atoms, temp_name):
    """
    Write an Orca input file. Currently, this function is redundant
    with the write_qchem function.

    A template file needs to be prepared which has XXX__POS__XXX
    instead of the atomic positions in the geometry.

    All other options must be specified in this file, e.g.
    symmetry, and the gradient keyword to get the gradients.

    Parameters
    ----------
    file_name : str
        Name of the Orca input file to be written
    atoms : list of Atom objects
        Atoms to be calculated with Orca
    temp_name : str
        Name of the template file
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")
    # replaces XXX__POS__XXX with atom coordinates
    for line in temp_content:
        if "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                out_file.write(atomStr)
        else:
            out_file.write(line)
    out_file.close()
    return

def write_mopac(file_name, atoms, temp_name):
    """
    Write a MOPAC input file.

    A template file needs to be prepared which has XXX__POS__XXX
    instead of the atomic positions below the keywords section.

    All other options also need to be specified in this file, e.g.
    

    Parameters
    ----------
    file_name : str
        Name of the QChem input file to be written
    atoms : list of Atom objects
        Atoms to be calculated with QChem
    temp_name : str
        Name of the template file
    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()

    out_file = open(file_name, "w")
    # replaces XXX__POS__XXX with atom coordinates and 1 after 
    # each coordinate to fulfill MOPAC requirements
    for line in temp_content:
        if "XXX__POS__XXX" in line:
            for atom in atoms:
                atomStr = "{:>6} {:10.6f} {:3} {:10.6f} {:3} {:10.6f} {:3}".format(
                    atom.elem, atom.x, 1, atom.y, 1, atom.z, 1) + "\n"
                out_file.write(atomStr)
        else:
            out_file.write(line)
    out_file.close()

    return  
  
def write_tinker_temp(file_name_1,file_name_2,reg1_atoms,reg2_atoms,point_charges=None):

    """
    Write a Tinker input template file to use it within a MOPAC calculation.

    No initial or template file is required in this case. However, Openbabel
    must be installed to use this function

    Parameters
    ----------
    file_name_1 : str
        Name of the Tinker .xyz input file to be written. This file contains
        the xyz coords of the QM region along with its connectivity, and the
        xyz coords of the QM'region
    file_name_2 : str
        Name of the Tinker .key input file to be written. This file contains 
        the semiempirical charges and the point charges and other parameters
    reg1_atoms : list of Atom objects
        Cartesian coords of the atoms in the QM region
    reg2_atoms : list of Atom objects
        Cartesian coords and charges of the atoms in the QM' region.
    point_charges : str
        Cartesian position and values of the high level point charges

    """
    if file_name_1=="rl.temp":
        with open(file_name_2) as temp_file:
            temp_content = temp_file.readlines()
        out_file_rl = open(file_name_1, "w") 
        for line in temp_content:
            if "XXX__FIX__XXX" in line:
                for atom in reg2_atoms:
                    atomStr = "{:>6} {:10.6f} 1 {:10.6f} 1 {:10.6f} 1 ".format(
                        atom.elem, atom.x, atom.y, atom.z) + "\n"
                    out_file_rl.write(atomStr)
            else:
                out_file_rl.write(line)
        out_file_rl.close()
    else:    
        # Delete the output files if they already exists
        # to avoid overwriting 
        if os.path.exists(file_name_1):
            subprocess.call("rm " + file_name_1, shell=True)
        if os.path.exists(file_name_2):
            subprocess.call("rm " + file_name_2, shell=True)

        # Create a .xyz file with the atoms of the QM region to get their
        # conectivity within the Tnk format using Openbabel
        write_xyz("init_geom.xyz",reg1_atoms)
   
        out_file = open("init_geom.xyz", "w")
        out_file.write(str(len(reg1_atoms)) + "\n") 
        out_file.write("" + "\n")
        for atom in reg1_atoms:
            atomStr = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                        atom.elem, atom.x, atom.y, atom.z) + "\n"
            out_file.write(atomStr)
        out_file.close()    

        try:
            subprocess.call("obabel -ixyz init_geom.xyz -otxyz -O init_geom.xyz", shell=True)
        except:
            pass

        with open("init_geom.xyz") as temp_file:
            temp_content = temp_file.readlines()

        out_file_xyz = open(file_name_1, "w")
        out_file_key = open(file_name_2, "w")

        n_total_atoms = str(len(reg1_atoms) + len(reg2_atoms))
        out_file_xyz.write(n_total_atoms + "\n")
        out_file_xyz.write("\n")

        for line in temp_content[1:]:
            out_file_xyz.write(line)

        cont = len(reg1_atoms) + 1
        for atom in reg2_atoms:
            atomStr = "{:>6} {:>5} {:10.6f} {:10.6f} {:10.6f} {:>5}".format(
                cont, "DM", atom.x, atom.y, atom.z, 42) + "\n"
            out_file_xyz.write(atomStr)
            cont += 1

        out_file_xyz.close()

        out_file_key.write("# Output Control" + "\n")
        out_file_key.write("DEBUG" + "\n")
        out_file_key.write("# Force Field Selection" + "\n")
        out_file_key.write("PARAMETERS    oplsaa.prm" + "\n")
        out_file_key.write("# Partial Structure" + "\n")
        out_file_key.write("GROUP    1 -1 " + str(len(reg1_atoms)) + "\n")
        out_file_key.write("GROUP    2 " + str(-(len(reg1_atoms)+1))+" "+str(n_total_atoms)+"\n")
        out_file_key.write("GROUP-SELECT    1 1 0.0" + "\n")
        out_file_key.write("GROUP-SELECT    1 2 1.0" + "\n")
        out_file_key.write("GROUP-SELECT    2 2 0.0" + "\n")
        out_file_key.write("" + "\n")

        cont = len(reg1_atoms) + 1
        if point_charges is not None:    
            for point in point_charges:
                point_str = "{:<8} {:<7} {:10.6f}".format(
                        "CHARGE", -cont, point.q) + "\n"
                out_file_key.write(point_str)
                cont += 1
        else:
            for atom in reg2_atoms:
                atomStr = "{:<8} {:<7} {:10.6f}".format(
                        "CHARGE", -cont, atom.q) + "\n"
                out_file_key.write(atomStr)
                cont += 1

        out_file_key.close()

        return

def write_tinker_xyz(file_name, atoms, temp_name):
    """
    Write a Tinker input file from a template to use within a MOPAC calculation.

    The template file is created automatically calling fro_prep_run.py. Openbabel 
    is required for this purpose. The atom types of the QM region must be adjusted 
    manually for the system under study

    Parameters
    ----------
    file_name : str
        Name of the Tinker input file to be written. This file contains
        the xyz coords of the QM region along with its connectivity, and the
        xyz coords of the QM'region
    atoms : list of Atom objects
        Cartesian coords of the atoms in the QM region
    temp_name : str
        Name of the template file        

    """
    with open(temp_name) as temp_file:
        temp_content = temp_file.readlines()
    out_file = open(file_name, "w")
    for line in temp_content[:2]:
        out_file.write(line)
    cont = 2
    for atom in atoms:
        temp_line = temp_content[cont].split()
        temp_line[1] = atom.elem
        temp_line[2] = str(atom.x)
        temp_line[3] = str(atom.y)
        temp_line[4] = str(atom.z)
        new_line = "   ".join(temp_line)
        out_file.write(new_line + "\n")
        cont += 1
    for line in temp_content[cont:]:
        out_file.write(line)
    out_file.close()
    return

def write_coord(in_atoms):
    """
    Write a Turbomole coord file

    The written units are in Bohr but we conserve Angstrom units in this program

    Parameters
    ----------
    in_atoms : Atom objects
        Atoms to write

    """
    bohrconv = 1.88973
    coord_file = open("coord", "w")
    coord_file.write("$coord\n")
    for atom in in_atoms:
        form_string = "{:10.6f} {:10.6f} {:10.6f} {:>6}\n".format(
            atom.x * bohrconv, atom.y * bohrconv, atom.z * bohrconv, atom.elem.lower())
        coord_file.write(form_string)
    coord_file.write("$end\n")
    return


def write_cube(in_name, origin, vectors, x_num, y_num, z_num, atoms, vals, comment=""):
    """
    Write a file in cube format

    Parameters
    ----------
    name : str
        Name of the output file
    origin : numpy array of length 3
        Origin if the cube grid
    vectors : 3x3 numpy array
        Vectors defining a voxel
    x_num, y_num, z_num : ints
        Numbers of voxels along each edge
    atoms : list of Atom objects
        Atoms to include in the cube file
    vals : numpy array
        The values to be entered at each voxel in the order x1, y1, z1, x1, y1, z2 etc.

    """
    bohrconv = 1.88973

    out_file = open(in_name, "w")

    # Header
    out_file.write("Cube file generated by fromage, units in Angstrom\n")
    out_file.write(comment + "\n")
    orig_line = "{:6d} {:10.6f} {:10.6f} {:10.6f}\n".format(
        len(atoms), origin[0] * bohrconv, origin[1] * bohrconv, origin[2] * bohrconv)
    out_file.write(orig_line)

    nums = [x_num, y_num, z_num]
    xyz_lines = ["{:6d} {:10.6f} {:10.6f} {:10.6f}\n".format(
        nums[i], vectors[i][0] * bohrconv, vectors[i][1] * bohrconv, vectors[i][2] * bohrconv) for i in list(range(3))]
    out_file.write(xyz_lines[0] + xyz_lines[1] + xyz_lines[2])

    # Atoms
    for atom in atoms:
        atom_line = "{:3d}{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n".format(
            atom.at_num, atom.at_num, atom.x * bohrconv, atom.y * bohrconv, atom.z * bohrconv)
        out_file.write(atom_line)

    # Values
    for i, val in enumerate(vals):
        out_file.write("{:10.6f}".format(val))
        if (i + 1) % 6 == 0:
            out_file.write("\n")

    return


def write_lat_vec(in_name, vectors):
    """
    Write vectors to a file

    Parameters
    ----------
    in_name : str
        Name of the file
    vectors : 3 x 3 numpy array
        Lattice vectors

    """
    out_file = open(in_name, "w")
    for line in [0, 1, 2]:
        out_file.write("{:15.11f}{:15.11f}{:15.11f}\n".format(
            vectors[line][0], vectors[line][1], vectors[line][2]))
    out_file.close()
