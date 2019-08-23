"""Defines the Calc objects which are used to run calculations

The purpose of the concrete classes is to handle all of the specificities of
each program so as to be able to call a calculation and read its output in
the same way from the main program regardless of which program it is.

As such if any clunky code is necessary due to an external program's specific
preferences, it should be contained here in its Calc object member methods.
"""
import subprocess
import os

from fromage.utils.mol import Mol
from fromage.io import edit_file as ef
from fromage.io import read_file as rf
from fromage.utils import handle_atoms as ha


bohrconv = 1.88973  # Something in Angstrom * bohrconv = Something in Bohr


def setup_calc(calc_name, calc_type):
    """
    Return a calculation of the correct subclass

    """
    calc_type = calc_type.lower()
    calc_types = {"gaussian" : Gauss_calc,
                  "gaussian_cas" : Gauss_CAS_calc,
                  "molcas" : Molcas_calc,
                  "turbomole" : Turbo_calc,
                  "turbomole_scf" : Turbo_SCF_calc,
                  "turbomole_tddft" : Turbo_calc_TDDFT,
                  "dftb" : DFTB_calc}
    try:
        out_calc = calc_types[calc_type](calc_name)
    except KeyError:
        print("Unercognised program: " + calc_type)

    return out_calc


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

    def read_out(self, positions, in_mol=None, in_shell=None):
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
        subdir = os.getcwd()
        os.chdir(self.here)
        with open("geom_mol.xyz", "a") as geom_m_file:
            geom_m_file.write(str(len(in_mol)) + "\n")
            geom_m_file.write(self.calc_name + "\n")
            for atom in ha.array2atom(in_mol, positions):
                atom_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_m_file.write(atom_str)
        # the inner and middle regions
        with open("geom_cluster.xyz", "a") as geom_c_file:
            geom_c_file.write(
                str(int((len(positions) / 3) + len(in_shell))) + "\n")
            geom_c_file.write(self.calc_name + "\n")
            for atom in ha.array2atom(in_mol, positions):
                atom_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_c_file.write(atom_str)
            for atom in in_shell:
                atom_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_c_file.write(atom_str)
        os.chdir(subdir)
        return


class DFTB_calc(Calc):
    """
    Calculation of DFTB+ tested with v18.2

    """

    def run(self, atoms):
        """
        Write a DFTB .gen file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with DFTB+
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method

        """
        dftb_path = os.path.join(self.here, self.calc_name)
        os.chdir(dftb_path)

        mol = Mol(atoms)

        region_2_file = "r2.xyz"
        if os.path.exists(region_2_file):
            mol_r2 = rf.mol_from_file(region_2_file)
            mol += mol_r2

        mol.write_xyz("geom.xyz")
        subprocess.call("xyz2gen geom.xyz", shell=True)
        # Run DFTB+
        proc = subprocess.Popen("dftb+ > dftb_out", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a DFTB+ detailed.out file while printing geometry updates

        To update the geom files, include in_mol and in_shell

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
        in_mol : list of Atom objects, optional
            Atoms in the inner region. Include to write geom files
        in_shell : list of Atom objects, optional
            Atoms in the middle region. Include to write geom files
        Returns
        -------
        energy : float
            Energy calculated by DFTB+ in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        dftb_path = os.path.join(self.here, self.calc_name)
        os.chdir(dftb_path)

        energy, gradients_b, scf_energy = rf.read_dftb_out("detailed.out")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)
        return (energy, gradients, scf_energy)


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

        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)

        ef.write_gauss(self.calc_name + ".com", atoms,
                       [], self.calc_name + ".temp")
        proc = subprocess.Popen(
            "${FRO_GAUSS} " + self.calc_name + ".com", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a Gaussian .chk file while printing geometry updates

        To update the geom files, include in_mol and in_shell

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
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
        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)

        # stdout=FNULL to not have to read the output of formchk
        FNULL = open(os.devnull, 'w')
        subprocess.call("formchk gck.chk", stdout=FNULL, shell=True)
        energy, gradients_b, scf_energy = rf.read_fchk("gck.fchk")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)

        return (energy, gradients, scf_energy)

    def read_out_mol(self, pop="EPS"):
        """Read the output log file and return Mol"""
        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)
        out_mol = rf.mol_from_gauss(self.calc_name + ".log")
        os.chdir(self.here)

        return out_mol

class Gauss_CAS_calc(Calc):
    """
    Calculation with Gaussian 09 for CAS calculations
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

        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)

        ef.write_gauss(self.calc_name + ".com", atoms,
                       [], self.calc_name + ".temp")
        proc = subprocess.Popen(
            "${FRO_GAUSS} " + self.calc_name + ".com", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a Gaussian .chk file while printing geometry updates

        To update the geom files, include in_mol and in_shell

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
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
        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)

        energy_e, grad_e, energy_g, grad_g = rf.read_g_cas(
            self.calc_name + ".log")
        # fix gradients units to Hartree/Angstrom
        grad_e = grad_e * bohrconv
        grad_g = grad_g * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        grad_e = grad_e[:len(positions)]
        grad_g = grad_g[:len(positions)]

        os.chdir(self.here)

        return (energy_e, grad_e, energy_g, grad_g)


def turbo_redefine(atoms):
    """Update Turbomole mos and run actual"""
    FNULL = open(os.devnull, 'w')
    ef.write_coord(atoms)
    # Update mos
    subprocess.call("rm -f mos", shell=True)
    with open("define_feed", "w") as tmp_def_in:
        # define input for Huckel guess
        tmp_def_in.write("\n\n\neht\n\n\n\n\n\n\n\n*\n\n")
    subprocess.call("define < define_feed", stdout=FNULL, shell=True)
    subprocess.call("rm -f define_feed", shell=True)
    subprocess.call("actual -r", shell=True)
    return

class Turbo_calc_TDDFT(Calc):
    """
    Calculation of TDDFT energy and gradients with Turbomole 7.0

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

        turbo_redefine(atoms)

        # Run Turbomole
        proc = subprocess.Popen(
            "dscf > dscf.out && egrad > grad.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a Turbomole grad.out file while printing geometry updates

        To update the geom files, include in_mol and in_shell

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
        in_mol : list of Atom objects, optional
            Atoms in the inner region. Include to write geom files
        in_shell : list of Atom objects, optional
            Atoms in the middle region. Include to write geom files
        Returns
        -------
        energy : float
            Energy calculated by Turbomole in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)

        energy, gradients_b, scf_energy = rf.read_tb_grout("grad.out")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)
        return (energy, gradients, scf_energy)


class Turbo_calc(Calc):
    """
    Calculation of CC2 with Turbomole 7.0

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

        turbo_redefine(atoms)

        # Run Turbomole
        proc = subprocess.Popen(
            "dscf > dscf.out && ricc2 > ricc2.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a Turbomole ricc2.out file while printing geometry updates

        To update the geom files, include in_mol and in_shell

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
        in_mol : list of Atom objects, optional
            Atoms in the inner region. Include to write geom files
        in_shell : list of Atom objects, optional
            Atoms in the middle region. Include to write geom files
        Returns
        -------
        energy : float
            Energy calculated by Turbomole in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)

        energy, gradients_b, scf_energy = rf.read_ricc2("ricc2.out")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)
        return (energy, gradients, scf_energy)


class Turbo_SCF_calc(Calc):
    """
    Calculation of SCF like DFT or HF with Turbomole 7.0

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

        turbo_redefine(atoms)

        # Run Turbomole
        proc = subprocess.Popen(
            "dscf > dscf.out && grad > grad.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a Turbomole gradient file while printing geometry updates

        To update the geom files, include in_mol and in_shell

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
        in_mol : list of Atom objects, optional
            Atoms in the inner region. Include to write geom files
        in_shell : list of Atom objects, optional
            Atoms in the middle region. Include to write geom files
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

        energy, gradients_b = rf.read_tbgrad("gradient")
        scf_energy = energy
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)
        return (energy, gradients, scf_energy)


class Molcas_calc(Calc):
    """
    Calculation with Molcas 8.0
    """

    def run(self, atoms):
        """
        Write a Molcas input file and return a subprocess.Popen

        Make sure the input file is called [name of calculation].input
        e.g. mh.input and the geometry file in Gateway is called geom.xyz

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Gaussian
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method

        """
        molcas_path = os.path.join(self.here, self.calc_name)
        os.chdir(molcas_path)

        # Write a temporary geom file for molcas to read
        ef.write_xyz("geom.xyz", atoms)

        proc = subprocess.Popen(
            "molcas molcas.input -f", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Analyse a Molcas .input file while printing geometry updates

        To update the geom files, include in_mol and in_shell. Also removes
        molcas.*

        Parameters
        ----------
        positions : list of floats
            List of atomic coordinates, important for truncation of gradients
            if too many are calculated
        in_mol : list of Atom objects, optional
            Atoms in the inner region. Include to write geom files
        in_shell : list of Atom objects, optional
            Atoms in the middle region. Include to write geom files
        Returns
        -------
        energy : float
            Energy calculated by Gaussian in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        molcas_path = os.path.join(self.here, self.calc_name)
        os.chdir(molcas_path)

        energy, gradients_b, scf_energy = rf.read_molcas("molcas.log")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)

        # for large molcas wavefunction information
        subprocess.call("rm -rf molcas.*", shell=True)

        return (energy, gradients, scf_energy)
