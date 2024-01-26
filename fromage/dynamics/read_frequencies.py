#!/usr/bin/env python
#################################################
#                  Author                       # 
#                                               #
#          Federico J. Hernandez                #
#                                               #
# Routine based on a similar one from PyRAI2MD  #
#            for gaussian and Orca              #
#                                               #
#################################################

import os,sys
import numpy as np
from numpy import linalg as la
from abc import ABC, abstractmethod
from fromage.dynamics.periodic_table import Element

class Normal_Modes():
    
    def __init__(self):
        self.n_frequencies = 0
        self.frequencies = []
        self.intensities = []
        self.n_atoms = 0
        self.atom_types = []
        self.atom_coords = []
        self.displacements = []
        self.reduced_mass = []
        self.atomic_mass = []
        self.atomic_charge = []


class File_Reader(ABC):

    @abstractmethod
    def read_file(self,input_file):
        pass


class Molden_Reader(File_Reader):

    def read_file(self, input_file):
        """
        """
        if os.path.isfile(input_file + ".freq.molden"):
            input_file = input_file + ".freq.molden"
        else:
            sys.exit("Input file {}.freq.molden does not exist".format(input_file))

        modes = Normal_Modes()

        with open(input_file, 'r') as rf:
            for n,line in enumerate(rf):
                if "[N_FREQ]" in line:
                    modes.n_frequencies = int(rf.readline())
                elif "[FREQ]" in line:
                    for i in range(modes.n_frequencies):
                        modes.frequencies.append(rf.readline().split())
                elif "[INT]" in line:
                    for i in range(modes.n_frequencies):
                        modes.intensities.append(rf.readline().split())
                elif "[NATOM]" in line:
                    modes.n_atoms = int(rf.readline())
                elif "[FR-COORD]" in line:
                    for i in range(modes.n_atoms):
                        modes.atom_coords.append(rf.readline().split())
                elif " vibration" in line:
                    vibration = []
                    for i in range(modes.n_atoms):
                        vibration.append(rf.readline().split())
                    modes.displacements.append(vibration)
                elif "[RMASS]" in line:
                    for i in range(modes.n_frequencies):
                        modes.reduced_mass.append(float(rf.readline()))

        for i,atom in enumerate(modes.atom_coords):
            modes.atom_types.append( atom[0] )
            modes.atomic_mass.append( Element(atom[0]).getMass() )
            modes.atomic_charge.append( Element(atom[0]).getNuc() )
            modes.atom_coords[i] = [float(atom[1]), float(atom[2]), float(atom[3])]

        modes.frequencies = np.array(modes.frequencies).astype(float).reshape(modes.n_frequencies,1)
        modes.intensities = np.array(modes.intensities).astype(float)
        modes.displacements = np.array(modes.displacements).astype(float).reshape(modes.n_frequencies, modes.n_atoms, 3)
        modes.atomic_mass = np.array(modes.atomic_mass).reshape(modes.n_atoms, 1)
        modes.atomic_charge = np.array(modes.atomic_charge)
        modes.atom_coords = np.array(modes.atom_coords)
        modes.reduced_mass = np.array(modes.reduced_mass).reshape(modes.n_frequencies, 1)

        for freq in modes.frequencies:
            if freq < 0:
                print("Imaginary frequency: {}".format(freq))
        modes.frequencies = np.absolute(modes.frequencies)

        return modes

        # write the eigenvectors to a file
        write_hess_eigenvecs(modes.displacements,modes.frequencies)


class Bagel_Reader(File_Reader):

    def read_file(self, input_file):
        pass

def count_atom(formula):
    na = ''
    natom = 0
    for x in formula + 'X':
        if x in '01234567890':
            na += x
        else:
            try:
                natom += int(na)
                na = ''
            except ValueError:
                pass

    return natom


def g16_xyz(cart_list):
    atoms = []
    xyz = []
    for line in cart_list:
        _, a, t, x, y, z = line.split()
        atoms.append(int(a))
        xyz.append([float(x), float(y), float(z)])
 
    return atoms, xyz


def g16_modes(vect_list):
    modes = []
    for block in vect_list:
        block = np.array([x.split()[2:] for x in block]).astype(float)
        nmode = int(block.shape[1] / 3)
        for n in range(nmode):
            modes.append(block[:, n: n + 3])

    modes = np.array(modes)

    return modes

def g16_format(data, dshape):
    ## formatting data
    dlist = []
    for i in data:
        dlist += [float(x) for x in i.split()]
    dlist = np.array(dlist)
    dlist = dlist.reshape(dshape)

    return dlist

class Gauss_Reader(File_Reader):
    
    def read_file(self, input_file):
        """
        Selects if it is going to read the normal modes
        eigenvectors from the log or the fchk file. 
        In the case of reading from the fchk file, the 
        option SaveNormalModes must be included in the
        Gaussian input file
        """
        if os.path.exists('%s.freq.fchk' % input_file):
            print('g16: read fchk and log')
            return self.read_g16_fchk(input_file)
        else:
            print('g16: read log only')
        return self.read_g16_log(input_file) 

    def read_g16_log(self,input_file):
        """
         This function reads .freq.g16 file (Gaussian .log file) and return all data 
         in the object modes
         The g16 saves normalized unmass-weighted normal modes
        """

        modes = Normal_Modes()
     
        with open('%s.freq.g16' % input_file, 'r') as raw:
            log = raw.read().splitlines()

        ## extracting data from log
        freqs_list = []
        rmass_list = []
        inten_list = []
        cart_list = []
        vect_list = []
        natom = 0
        reading=False
        for n, line in enumerate(log):

            if reading:
                if len(line.split()) == 4:
                    natom += 1
                else:
                    reading = False

            if "Charge" in line and "Multiplicity" in line:
                reading = True

            if 'Input orientation' in line:
                cart_list = log[n + 5: n + 5 + natom]

            if 'Frequencies -- ' in line:
                f = line.split(' -- ')[-1]
                freqs_list.append(f)

            if 'Red. masses -- ' in line:
                r = line.split(' -- ')[-1]
                rmass_list.append(r)

            if 'IR Inten    --' in line:
                i = line.split(' -- ')[-1]
                inten_list.append(i)

            if ' Atom  AN ' in line:
                v = log[n + 1: n + 1 + natom]
                vect_list.append(v)

        atoms, xyz = g16_xyz(cart_list)
        atoms = np.array([Element(str(int(i))).getSymbol() for i in atoms])

        freqs = g16_format(freqs_list, [-1, 1])
        rmass = g16_format(rmass_list, [-1, 1])
        inten = g16_format(inten_list, [-1])

        nmodes = g16_modes(vect_list)
        nmodes2print = nmodes.copy() # copy the normalised eigenvectors (non mass weighted)

        amass = [Element(i).getMass() for i in atoms]
        amass = np.array(amass)
        amass = amass.reshape((natom, 1))

        achrg = [Element(i).getNuc() for i in atoms]
        achrg = np.array(achrg)
        achrg = achrg.reshape((natom, 1))

        nmodes = np.array([i / la.norm(i * amass ** 0.5) for i in nmodes]) # convert to normalized 

        modes.n_atoms = natom
        modes.atom_types = atoms
        modes.atomic_mass = amass
        modes.reduced_mass = rmass
        modes.atomic_charge = achrg
        modes.atom_coords = np.array(xyz) 
        modes.n_frequencies = len(freqs)
        modes.frequencies = freqs 
        modes.intensities = inten
        modes.displacements = nmodes

        for freq in modes.frequencies:
            if freq < 0:
                print("Imaginary frequency: {}".format(freq))
        modes.frequencies = np.absolute(modes.frequencies)

        # write the eigenvectors to a file
        write_hess_eigenvecs(nmodes2print,modes.frequencies)

        return modes

    def read_g16_fchk(self,input_file):
        """
         This function reads .freq.g16 (Gaussian .log file) and the .freq.fchk files and return
         all data in the object modes
         The g16 saves normalized unmass-weighted normal modes
        """

        modes = Normal_Modes()

        with open('%s.freq.g16' % input_file, 'r') as raw:
            log = raw.read().splitlines()
        with open('%s.freq.fchk' % input_file, 'r') as raw:
            fchk = raw.read().splitlines()

        ## extracting data from log
        freqs_list = []
        rmass_list = []
        inten_list = []

        for line in log:
            if 'Frequencies -- ' in line:
                f = line.split(' -- ')[-1]
                freqs_list.append(f)

            if 'Red. masses -- ' in line:
                r = line.split(' -- ')[-1]
                rmass_list.append(r)

            if 'IR Inten    --' in line:
                i = line.split(' -- ')[-1]
                inten_list.append(i)

        freqs = g16_format(freqs_list, [-1, 1])
        rmass = g16_format(rmass_list, [-1, 1])
        inten = g16_format(inten_list, [-1])

    ## extracting data from fchk

        natom = 0
        atom_list = []
        cart_list = []
        nmode = 0
        vect_list = []

        for n, line in enumerate(fchk):
            if 'Atomic numbers' in line:
                natom = int(line.split()[-1])
                atom_line = int(natom / 6) + (natom % 6 > 0)
                atom_list = fchk[n + 1: n + 1 + atom_line]

            if 'Current cartesian coordinates' in line:
                ncart = int(line.split()[-1])
                cart_line = int(ncart / 5) + (ncart % 5 > 0)
                cart_list = fchk[n + 1: n + 1 + cart_line]

            if 'Number of Normal Modes' in line:
                nmode = int(line.split()[-1])

            if 'Vib-Modes' in line:
                nvect = int(line.split()[-1])
                vect_line = int(nvect / 5) + (nvect % 5 > 0)
                vect_list = fchk[n + 1: n + 1 + vect_line]

        atoms = g16_format(atom_list, [-1])
        atoms = np.array([Element(str(int(i))).getSymbol() for i in atoms])
 
        xyz = g16_format(cart_list, [natom, 3])
        nmodes = g16_format(vect_list, [nmode, natom, 3])
        nmodes2print = nmodes.copy() # copy the normalised eigenvectors (non mass weighted)

        amass = [Element(i).getMass() for i in atoms]
        amass = np.array(amass)
        amass = amass.reshape((natom, 1))

        achrg = [Element(i).getNuc() for i in atoms]
        achrg = np.array(achrg)
        achrg = achrg.reshape((natom, 1))

        nmodes = np.array([i / la.norm(i * amass ** 0.5) for i in nmodes])  # convert to normalized 
                                                                            # mass-weighted

        modes.n_atoms = natom
        modes.atom_types = atoms #[Element(str(int(i))).getSymbol() for i in atoms]
        modes.atomic_mass = amass
        modes.reduced_mass = rmass
        modes.atomic_charge = achrg
        modes.atom_coords = np.array(xyz)
        modes.n_frequencies = len(freqs)
        modes.frequencies = freqs
        modes.intensities = inten
        modes.displacements = nmodes

        for freq in modes.frequencies:
            if freq < 0:
                print("Imaginary frequency: {}".format(freq))
        modes.frequencies = np.absolute(modes.frequencies)

        # write the eigenvectors to a file
        write_hess_eigenvecs(nmodes2print,modes.frequencies)

        return modes

class Orca_Reader(File_Reader):

    def read_file(self, input_file):
        """
        Function to read the .hess file from Orca and return all 
        data in the modes object. Orca saves normalised unmass-weighted 
        normal modes
        """

        modes = Normal_Modes()

        with open('%s.freq.orca' % input_file, 'r') as orca:
            hess = orca.read().splitlines()

        atoms = np.zeros(0)
        xyz = np.zeros(0)
        freqs = np.zeros(0)
        inten = np.zeros(0)
        nmodes = np.zeros(0)
 
        ## extrac data from .hess file

        for n, line in enumerate(hess):
            if '$vibrational_frequencies' in line:
                natom = int(int(hess[n + 1]) / 3)
                f = hess[n + 2: n + 2 + natom * 3]
                freqs = np.array([x.split() for x in f]).astype(float)[:, 1].reshape((-1, 1))
            if '$ir_spectrum' in line:
                natom = int(int(hess[n + 1]) / 3)
                i = hess[n + 2:n + 2 + natom * 3]
                inten = np.array([x.split() for x in i]).astype(float)[:, 1].reshape((-1, 1))
            if '$atoms' in line:
                natom = int(hess[n + 1])
                coord = hess[n + 2: n + 2 + natom]
                coord = np.array([x.split() for x in coord])
                atoms = coord[:, 0]
                xyz = coord[:, 2: 5].astype(float)
            if '$normal_modes' in line:
                nmode = int(hess[n + 1].split()[0])
                nline = (nmode + 1) * (int(nmode / 5) + (nmode % 5 > 0))
                vects = hess[n + 2: n + 2 + nline]
                nmodes = [[] for _ in range(nmode)]
                for m, i in enumerate(vects):
                    row = m % (nmode + 1) - 1
                    if row >= 0:
                        nmodes[row] += [float(j) for j in i.split()[1:]]
                nmodes = np.array(nmodes).T.reshape((nmode, int(nmode / 3), 3))  # Transpose array


        # filter trans-rot freqs and modes
        realfreq = []
        for n, i in enumerate(freqs):
            if i > 0:
                realfreq.append(n)
            elif i < 0:
                print('imaginary frequency: %10f mode: %6d (converted to real)' % (i, n + 1))  
                realfreq.append(n)

        nmode = len(realfreq)
        amass = [Element(atm).getMass() for atm in atoms]
        amass = np.array(amass)
        amass = amass.reshape((natom, 1))

        achrg = [Element(atm).getNuc() for atm in atoms]
        achrg = np.array(achrg)
        achrg = achrg.reshape((natom, 1))

        modes.n_atoms = len(atoms)
        modes.atom_types = atoms #[Element(str(int(i))).getSymbol() for i in atoms]
        modes.atomic_mass = amass
        modes.reduced_mass = np.array([0 for _ in range(nmode)]).reshape((-1, 1))
        modes.atomic_charge = achrg
        modes.atom_coords = np.array(xyz)
        modes.frequencies = freqs[realfreq].reshape((-1, 1))
        modes.n_frequencies = len(modes.frequencies)
        modes.intensities = inten[realfreq]
        modes.displacements = nmodes[realfreq]

        for freq in modes.frequencies:
            if freq < 0:
                print("Imaginary frequency: {}".format(freq))
        modes.frequencies = np.absolute(modes.frequencies)

        # write the eigenvectors to a file
        #write_hess_eigenvecs(nmodes,modes.frequencies)

        return modes

class Turbomole_Reader(File_Reader):

    def read_file(self, input_file):
        """
        """
        bohr_rad = 0.52917720859
        reading = False
        reading_mode = False
        atom_coords = []
        frequencies = []
        intensities = []
        reduced_mass = []
        x_mode = []
        y_mode = []
        z_mode = []
        displacements = []
       
        with open('%s.freq.turbomole' % input_file, 'r') as data:
            turbo = data.read().splitlines()

            ## extracting data from TM output
        for n,line in enumerate(turbo):
            if "atomic coordinates" in line:
                reading = True
                continue
#            if "center of nuclear mass" in line:
#                n_atoms = int(len(atom_coords))
#                continue
            if reading:
                if line.strip():
                    atom_coords.append(line.split()[0:4])
                else:
                    n_atoms = int(len(atom_coords))
                    reading = False
            if line.strip():
                if line.split()[0] == "frequency":
                    nums = [float(i) for i in line.split()[1:]]
                    frequencies.extend(nums)
                if line.split()[0] == "intensity" and line.split()[1] == "(km/mol)":
                    nums = [float(i) for i in line.split()[2:]]
                    intensities.extend(nums)
                if line.split()[0] == "RAMAN":
                    reading_mode = True
                    continue
                if line.split()[0] == "reduced" and line.split()[1] == "mass(g/mol)":
                    reading_mode = False
                    nums = [float(i) for i in line.split()[2:]]
                    reduced_mass.extend(nums)
                if reading_mode:
                    if line.strip():
                        if line.split()[2] =="x":
                            nums = [float(i) for i in line.split()[3:]]
                            x_mode.extend(nums)
                        if line.split()[0] =="y":
                            nums = [float(i) for i in line.split()[1:]]
                            y_mode.extend(nums)
                        if line.split()[0] =="z":
                            nums = [float(i) for i in line.split()[1:]]
                            z_mode.extend(nums)
        
        for x, y, z in zip(x_mode, y_mode, z_mode): # Asi guarda los x y z de los 6 primeros modos para el atomo 1, luego de los 6
                                            # primeros para el atomo 2, etc.CORREGIR!!
            displacements.append([x,y,z])
        mode_to_delete = []
        for i in range(len(frequencies)):
            if frequencies[i] == 0.:
                mode_to_delete.append(i)
            if frequencies[i]:
                displacements_2 = displacements
        for mode in sorted(mode_to_delete, reverse = True ):
            del frequencies[mode]
            del intensities[mode]
            del reduced_mass[mode]
            ini_atom = n_atoms * mode
            for i in sorted(range(n_atoms),reverse=True):
                del displacements[ini_atom + i]
        n_modes = int(3. * n_atoms - 6.)
        displacements = np.array(displacements).reshape(n_modes,n_atoms,3)
 
        print("n_modes")
        print(n_modes)
        print("")
        print(displacements[0,:,:])
        print("")
        print("")
        print("")
        print(displacements[-1,:,:])


        return 
#        pass
        
def write_hess_eigenvecs(nmodes,freqs):
    """
     This function writes the mass-weigthed hessian eigenvectors to a 
      eigenvectors.dat file and the normal modes frequencies to  
      frequencies.dat file

    """
    # Assuming nmodes is your original numpy array with the shape (nmodes, natoms, 3)
    # Replace this with the actual code that generates or reads nmodes

    # Transforming nmodes to new_array with shape (3*natoms, nmodes)
    new_array = nmodes.transpose(1, 2, 0).reshape(-1, nmodes.shape[0])

    header_comment = "Hessian mass-weighted eigenvectors ordered as (3*{natoms}, {nmodes_minus_6})"
    formatted_comment = header_comment.format(natoms=nmodes.shape[1], nmodes_minus_6=nmodes.shape[0] - 6)

    np.savetxt("eigenvectors.dat", new_array, header=formatted_comment, comments='# ')

    header_comment = "Vibrational Frequencies in cm-1 ({nmodes} Freqs)"
    formatted_comment = header_comment.format(nmodes=freqs.shape[0])

    np.savetxt("frequencies.dat", freqs, header=formatted_comment, comments='# ')

    return None 
