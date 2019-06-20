""" Read and extracts information from gaussian 09 log files
"""
from fromage.utils.exci_coupling import elements
import numpy as np
au2ev=27.211396132
def read_xyz(g09_file):
    """
    Opens a g09 log file and returns the first geometry in Input orientation.

    Iterators are used so that the file is not all loaded into memory, which
    can be expensive.

    The function searches for the following text pattern in the log file:

>                            Input orientation:
>    ---------------------------------------------------------------------
>    Center     Atomic      Atomic             Coordinates (Angstroms)
>    Number     Number       Type             X           Y           Z
>    ---------------------------------------------------------------------

    And will save the coordinates and the atomic symbols succeeding it

    Parameters
    ----------
    g09_file: Path to g09 log file
        File path

    Returns
    -------
    coordinates: List of lists
        Outer list is whole xyz file, each inner list is a line of the file containing
        the symbol and x,y,z coordinates

    """
    with open(g09_file) as f:
        # Get the number of atoms so we can iterate without loading the file into memory
        for line in f:
            # Ensures line is not blank
            if line.strip():
                if line.split()[0]=="NAtoms=":
                    natoms=(int(line.split()[3]))
                    break
        # Will hold the coordinates and symbols
        coordinates=[]
        # Reset the iterator to the top of the file
        f.seek(0)
        for line in f:
            if line.strip():
                if "Input orientation:" in line:
                    for i in range(5):
                        # Skip 5 lines to start of coordinates
                        line=next(f)
                    for i in range(natoms):
                        linesplit=line.split()

                        symb=str(elements.element(int(linesplit[1])).symbol)
                        x=float(linesplit[3])
                        y=float(linesplit[4])
                        z=float(linesplit[5])
                        coordinates.append([symb,x,y,z])
                        line=next(f)
                    break
                    f.close()
        return coordinates

def read_natoms(g09_file):
    """
    Opens a g09 log file and returns the number of atoms in the system

    Parameters
    ----------
    g09_file: Path to g09 log file
        File path

    Returns
    -------
    natoms: Integer

    """
    with open(g09_file) as f:
        # Get the number of atoms so we can iterate without loading the file into memory
        for line in f:
            # Ensures line is not blank
            if line.strip():
                if line.split()[0]=="NAtoms=":
                    natoms=(int(line.split()[3]))
                    break
        return natoms

def read_TD(g09_file,state):
    """
    Reads a G09 logfile and returns the Transition Dipole vector for the specified
    electronic state

    Parameters
    ----------
    g09_file: Path to g09 log file
        File path

    Returns
    -------
    TD: np.array
        1x3 array of x,y,z components of TD vector

    """
    with open(g09_file) as f:
        # Get the number of atoms so we can iterate without loading the file into memory
        for line in f:
            # Ensures line is not blank
            if line.strip():
                if " Ground to excited state transition electric dipole moments (Au):" in line:
                    for i in range(state+1):
                        line=next(f)
                    s,X,Y,Z,DipS,Osc = line.split()
                    break
        f.close()
        return np.array([float(X),float(Y),float(Z)])

def read_NTO(g09_file,natoms):
    """
    Reads a G09 logfile and returns the atomic centred Natural Transition Charges, obtained
    via the G09 input line:
    td=(nstates=1) nosymm Pop=NTO Density=(Transition=1)

    Parameters
    ----------
    g09_file: Path to g09 log file
        File path
    natoms: Integer
        Number of atoms
    Returns
    -------
    NTO: np.array
        N array of NTO charges in order of atomic positions

    """
    NTO=np.zeros(natoms)
    with open(g09_file) as f:
        # Get the number of atoms so we can iterate without loading the file into memory
        for line in f:
            # Ensures line is not blank
            if line.strip():
                if " Mulliken charges:" in line:
                    line = next(f)
                    line = next(f)
                    for i in range(natoms):
                        charge_line = line.split()
                        symbol = charge_line[1]
                        charge = float(charge_line[2])
                        # NTO charge is atomic number - charge
                        NTO[i] = elements.element(symbol).atomic-float(charge)
                        line = next(f)

        f.close()
    return NTO

def read_SCF(g09_file):
    """
    Opens a g09 log file and returns the final SCF energy

    Parameters
    ----------
    g09_file: Path to g09 log file
        File path

    Returns
    ----------
    SCF: float
        Final SCF energy

    """
    energies=[]
    with open(g09_file) as f:
        for line in f:
            # Ensures line is not blank
            if line.strip():
                if " SCF Done:" in line:
                    energies.append(float(line.split()[4]))
        f.close()
    SCF=energies[-1]
    return SCF

def read_ES(g09_file,state):
    """
    Opens a g09 log file and returns the energy difference between the ground and
    specified state in atomic units

    Parameters
    ----------
    g09_file: Path to g09 log file
        File path
    state: Integer
        Excited state number (<1)
    Returns
    -------
    ES: float
        Energy difference in atomic units between the ground and specified state

    """
    if state <1 :
        exit("Specified state must be an excited state (>0)")
    SCF=read_SCF(g09_file)
    energies=[]
    with open(g09_file) as f:
        for line in f:
            # Ensures line is not blank
            if line.strip():
                if " Excited State   {}".format(state) in line:
                    energies.append(float(line.split()[4]))
        f.close()

        ES = energies[-1]/au2ev

    return ES
