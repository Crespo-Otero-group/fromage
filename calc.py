"""Defines the Calc objects which are used to run calculations

The purpose of the concrete classes is to handle all of the specificities of
each program so as to be able to call a calculation and read its output in
the same way from the main program regardless of which program it is.

As such if any clunky code is necessary due to an external program's specific
preferences, it should be contained here in its Calc object member methods.
"""

import edit_file as ef
import read_file as rf
import handle_atoms as ha
import subprocess
import os

bohrconv = 1.88973  # Something in Angstrom * bohrconv = Something in Bohr


class Calc(object):
    """
    Abstract class for calculation objects

    Attributes
    ----------
    calc_name : str
        Name of the calculation, typically rl, ml, mh or mg
    """

    def __init__(self, calc_name_in=None, in_here=os.getcwd()):
        """Constructor which sets the calculation name"""
        self.calc_name = calc_name_in
        self.here = in_here

    def run(self, atoms):
        """
        Write all of the variable inputs necessary for one calculations

        """
        raise NotImplementedError("Please Implement this method")

    def read_out(self, positions, print_bool, in_mol, in_shell):
        """
        Read the output of the calculation and sometimes updates the geom_*.xyz files

        """
        raise NotImplementedError("Please Implement this method")

    def update_geom(self, positions, in_mol, in_shell):
        """
        Update the geom_mol.xyz and geom_cluster.xyz files

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates
        in_mol : list of Atom objects
            Atoms in the inner region
        in_shell : list of Atom objects
            Atoms in the middle region
        """
        with open("geom_mol.xyz", "a") as geom_m_file:
            geom_m_file.write(str(len(in_mol)) + "\n")
            geom_m_file.write(self.calc_name + "\n")
            for atom in ha.array2atom(in_mol, positions):
                atom_str = "{:>6} {:10.9f} {:10.6f} {:10.9f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_m_file.write(atom_str)
        # the inner and middle regions
        with open("geom_cluster.xyz", "a") as geom_c_file:
            geom_c_file.write(
                str((len(positions) / 3) + len(in_shell)) + "\n")
            geom_c_file.write(self.calc_name + "\n")
            for atom in ha.array2atom(in_mol, positions):
                atom_str = "{:>6} {:10.9f} {:10.6f} {:10.9f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_c_file.write(atom_str)
            for atom in in_shell:
                atom_str = "{:>6} {:10.9f} {:10.9f} {:10.9f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_c_file.write(atom_str)
        return


class Gauss_calc(Calc):
    """
    Calculation with Gaussian 09
    """

    def run(self, atoms):
        """
        Write a Gaussian input file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Gaussian
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method
        """
        ef.write_gauss(self.calc_name, self.calc_name + ".com",
                       atoms, [], self.calc_name + ".temp")
        proc = subprocess.Popen("g09 " + self.calc_name + ".com", shell=True)

        return proc

    def read_out(self, print_bool, positions=None, in_mol=None, in_shell=None):
        """
        Analyse a Gaussian .chk file while printing geometry updates

        Parameters
        ----------
        print_bool : bool
            If true, update the geometry files with the one found in .chk
        positions : list of floats, optional
            List of atomic coordinates
        in_mol : list of Atom objects, optional
            Atoms in the inner region
        in_shell : list of Atom objects, optional
            Atoms in the middle region
        Returns
        -------
        energy : float
            Energy calculated by Gaussian in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        # stdout=FNULL to not have to read the output of formchk
        FNULL = open(os.devnull, 'w')
        subprocess.call("formchk " + self.calc_name +
                        ".chk", stdout=FNULL, shell=True)
        energy, gradients_b, scf_energy = rf.read_fchk(
            self.calc_name + ".fchk")
        # fix gradients units
        gradients = gradients_b * bohrconv
        # update the geometry log
        if print_bool == True:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        return (energy, gradients, scf_energy)


class Turbo_calc(Calc):
    """
    Calculation with Turbomole 7.0
    """

    def run(self, atoms):
        """
        Write a Turbomole coord file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Gaussian
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method

        """
        FNULL = open(os.devnull, 'w')

        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)

        ef.write_coord(atoms)
        proc = subprocess.Popen(
            "dscf > dscf.out && ricc2 > ricc2.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, print_bool, positions=None, in_mol=None, in_shell=None):
        """
        Analyse a Turbomole ricc2.out file while printing geometry updates

        Parameters
        ----------
        print_bool : bool
            If true, update the geometry files with the one found in .chk
        positions : list of floats, optional
            List of atomic coordinates
        in_mol : list of Atom objects, optional
            Atoms in the inner region
        in_shell : list of Atom objects, optional
            Atoms in the middle region
        Returns
        -------
        energy : float
            Energy calculated by Gaussian in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)

        energy, gradients_b, scf_energy = rf.read_ricc2("ricc2.out")
        # fix gradients units
        gradients = gradients_b * bohrconv
        # update the geometry log
        if print_bool == True:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)
        return (energy, gradients, scf_energy)
