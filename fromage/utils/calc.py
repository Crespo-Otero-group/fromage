"""Defines the Calc objects which are used to run calculations

The purpose of the concrete classes is to handle all of the specificities of
each program so as to be able to call a calculation and read its output in
the same way from the main program regardless of which program it is.

As such if any clunky code is necessary due to an external program's specific
preferences, it should be contained here in its Calc object member methods.
"""
import numpy as np
import subprocess
import os

from fromage.utils.mol import Mol
from fromage.io import edit_file as ef
from fromage.io import read_file as rf
from fromage.utils import array_operations as ao


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
                  "turbomole_mp2" : Turbo_calc_MP2, 
                  "turbomole_scf" : Turbo_SCF_calc,
                  "turbomole_tddft" : Turbo_calc_TDDFT,
                  "dftb" : DFTB_calc,
                  "xtb" : xtb_calc,
                  "xtb_gfnff" : xtb_calc_gfnff,
                  "qchem" : Qchem,
                  "nwchem_dft": nwchem_calc_DFT,
                  "mopac": fomo_ci_calc,
                  "fomo-ci": fomo_ci_calc,
                  "orca": Orca_calc}
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

    def run(self, atoms, nprocs):
        """
        Write all of the variable inputs necessary for one calculations

        """
        raise NotImplementedError("Please Implement this method")
  
    def run_freq(self, atoms, nprocs):
        """
        Write all of the variable inputs necessary for a frequency calculation

        """
        raise NotImplementedError("Please Implement this method")

    def read_out(self, positions, in_mol=None, in_shell=None):
        """
        Read the output of the calculation and sometimes updates the geom_*.xyz files

        """
        raise NotImplementedError("Please Implement this method")

    def read_hessian(self):
        """
        Read the Hessian matrix elements
        """
        raise NotImplementedError("Please Implement this method")

    def read_nacs(self):
        """
        Read the Nonadiabatic coupling matrix elements    

        """
        raise NotImplementedError("Please Implement this method")

    def read_mu(self):
        """
        Read dipole vectors and dipole derivative matrix
        """
        raise NotImplementedError("Please Implement this method")
 
    def save_checkpoint(self):
        """
        Save all the importat info for a dynamics every "chk_stp" steps         

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
            for atom in ao.array2atom(in_mol, positions):
                atom_str = "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(
                    atom.elem, atom.x, atom.y, atom.z) + "\n"
                geom_m_file.write(atom_str)
        # the inner and middle regions
        with open("geom_cluster.xyz", "a") as geom_c_file:
            geom_c_file.write(
                str(int((len(positions) / 3) + len(in_shell))) + "\n")
            geom_c_file.write(self.calc_name + "\n")
            for atom in ao.array2atom(in_mol, positions):
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
    Calculation of DFTB+ tested with v22.2

    """
    def run(self, atoms, points_flex = None, nprocs=None):
        """
        Runs a DFTB+ force calculation using a .xyz file 
        and return a subprocess.Popen

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

        ef.write_dftb("geom.xyz",atoms,
                            [], self.calc_name + ".temp")
        if points_flex is not None:
            ef.write_dftb_charges("charges.dat", points_flex)
        # Run DFTB+
        proc = subprocess.Popen("dftb+ > dftb_out", shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, points_flex = None, nprocs=None):
        """
        Runs a DFTB+ freq calculation using a .xyz file 
        and return a subprocess.Popen

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

        ef.write_dftb("geom.xyz",atoms,
                            [], self.calc_name + ".temp", freq = True)
        if points_flex is not None:
            ef.write_dftb_charges("charges.dat", point_flex)
        # Run DFTB+
        proc = subprocess.Popen("dftb+ > dftb_out", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
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

        energy, gradients_bohr, scf_energy = rf.read_dftb_out("detailed.out")

        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long and fix gradients units to Hartree/Angstrom
        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                gradients = np.zeros(dim_flex)
            else:
                gradients = np.zeros(len(positions))
            # Fix gradients units to Hartree/Angstrom
            gradients[:len(positions)] = gradients_bohr[:len(positions)] * bohrconv
        else:
            # Fix gradients units to Hartree/Angstrom
            gradients = gradients_bohr[:len(positions)] * bohrconv

        os.chdir(self.here)
        return (energy, gradients, scf_energy)

    def read_charges(self):
        """
        Get the atomic charges of the whole system
        
        Returns
        ----------
        charges : array of atom charges
        """
        ## Add section to read CM5 charges
   
        dftb_path = os.path.join(self.here, self.calc_name)
        os.chdir(dftb_path)
        charges = rf.read_dftb_charges("detailed.out")
        os.chdir(self.here)
        return charges

    def read_mu(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Read dipole moment and dipole derivatives from a dftb+ output

        Returns
        ----------
        d_mu : array of dipole derivatives
        """

        dftb_path = os.path.join(self.here, self.calc_name)
        os.chdir(dftb_path)
        d_mu_tmp = rf.read_dftb_mu('born.out')

        #truncate the dipole derivatives matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3 * natoms_flex)
                d_mu = np.zeros((dim_flex,3))
            else:
                d_mu = np.zeros((len(positions),3))
            d_mu[:len(positions),:3] = d_mu_tmp[:len(positions),:3]
        else:
            d_mu = d_mu_tmp[:len(positions),:3]
        
        os.chdir(self.here)

        return d_mu

    def read_hessian(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Get the Hessian matrix from a dftb+ output

        Returns
        ----------
        hess : 3Natoms x 3Natoms array where Natoms is the amount of atoms in the
        QM region or plus the atoms in the flexible QM' region.
        """
        dftb_path = os.path.join(self.here, self.calc_name)
        os.chdir(dftb_path)
        hess_tmp = rf.read_hessian_dftb("hessian.out")
        #truncate the hessian matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                hess = np.zeros((dim_flex,dim_flex))
            else:
                hessian = np.zeros((len(positions),len(positions)))
            # Fix gradients units to Hartree/Angstrom
            hess[:len(positions),:len(positions)] = hess_tmp[:len(positions),:len(positions)]
        else:
            hess = hess_tmp[:len(positions),:len(positions)]

        os.chdir(self.here)
        return hess

class Gauss_calc(Calc):
    """
    Calculation with Gaussian 09/16
    """
    def run(self, atoms, point_flex = None, nprocs = None, state=None, states=None, singlestate=0, nac_coupling=[], soc_coupling=[]):
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

        if point_flex is not None:
            if state != None and states != None:
                ef.write_gauss(self.calc_name + ".com", atoms,
                           point_flex, self.calc_name + ".temp", state = state, states = states)
            else:
                ef.write_gauss(self.calc_name + ".com", atoms,
                           point_flex, self.calc_name + ".temp")
        else:
            if state != None and states != None:
                ef.write_gauss(self.calc_name + ".com", atoms,
                       [], self.calc_name + ".temp", state = state, states = states)
            else:
                ef.write_gauss(self.calc_name + ".com", atoms,
                       [], self.calc_name + ".temp")
        proc = subprocess.Popen(
            "${FRO_GAUSS} " + self.calc_name + ".com", shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, point_flex = None, nprocs=None):
        """
        Write a Gaussian input file for a normal modes calculation 
        and return a subprocess.Popen

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

        freq = True
        if point_flex is not None:
            ef.write_gauss(self.calc_name + ".com", atoms,
                       point_flex, self.calc_name + ".temp", freq = freq)
        else:
            ef.write_gauss(self.calc_name + ".com", atoms,
                       [], self.calc_name + ".temp", freq = freq)
        proc = subprocess.Popen(
            "${FRO_GAUSS} " + self.calc_name + ".com", shell=True)

        os.chdir(self.here)

        return proc


    def read_out(self, 
                 positions, 
                 dyn_bool=False, 
                 in_mol=None, 
                 in_shell=None, 
                 natoms_flex=None,
                 natoms = None,
                 state = None,
                 states = None,
                 mult = [],
                 singlestate = 0,
                 soc_coupling = []):
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
        proc_fchk = subprocess.call("formchk -0 gck.chk gck.fchk", shell=True)
        fchk_file="gck.fchk"

        nac = []
        soc = []

        if state is not None and states is not None:
            energy, gradients_b, scf_energy, nac, soc = rf.read_gauss_dyn(self.calc_name+".log",
                                                                          fchk_file,
                                                                          natoms,
                                                                          state,
                                                                          states,
                                                                          mult,
                                                                          singlestate,
                                                                          soc_coupling)

            # fix gradients units to Hartree/Angstrom
            gradients = gradients_b * bohrconv
#                                                                     #
#           ADD THE LINES TO ACCOUNT FOR THE FELIXIBILITY OF REGION 2 #
#                                                                     #


#        if dyn_bool:
#            #HERE I NEED TO CREATE A FUNCTION TO READ THE GRADIENTS AND THE ENERGIES
#            #SPECIFICALLY FOR THE DYNAMICS
#            energy = []
#            gradients_b = np.array([])
#            proc_fchk = subprocess.call("formchk -0 gck_GS.chk gck_GS.fchk", shell=True)
#            proc_fchk = subprocess.call("formchk -0 gck_ES.chk gck_ES.fchk", shell=True)
#            energy_ES, gradients_b_ES, scf_energy = rf.read_fchk("gck_ES.fchk")
#            energy_GS, gradients_b_GS, scf_energy = rf.read_fchk("gck_GS.fchk")
#            energy.append(scf_energy)
#            energy.append(energy_ES)
#            energy = np.reshape(np.array(energy),(-1,1))
#            gradients_b = np.concatenate((gradients_b_GS, gradients_b_ES))
#            # fix gradients units to Hartree/Angstrom
#            gradients = gradients_b * bohrconv        

        else:
            # stdout=FNULL to not have to read the output of formchk
            # FNULL = open(os.devnull, 'w')
#            proc_fchk = subprocess.call("formchk -0 gck.chk gck.fchk", shell=True)
            energy, gradients_b, scf_energy = rf.read_fchk("gck.fchk")
            # fix gradients units to Hartree/Angstrom
#            gradients = gradients_b * bohrconv
            # update the geometry log
            if in_mol != None:
                self.update_geom(positions, in_mol, in_shell)                      
            # truncate gradients if too long
#          
            # truncate gradients if too long and fix gradients units to Hartree/Angstrom
            if natoms_flex is not None:
                if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                    dim_flex = int(len(positions) + 3. * natoms_flex)
                    gradients = np.zeros(dim_flex)
                else:
                    gradients = np.zeros(len(positions))
                gradients[:len(positions)] = gradients_b[:len(positions)] * bohrconv
            else:
                gradients = gradients_b[:len(positions)] * bohrconv
##############################3 OLD IDEA ############################################## 
#            dim_flex = int(len(positions) + 3. * natoms_flex)
#            gradients = np.zeros(dim_flex)
#            gradients[:len(positions)] = gradients_b[:len(positions)] * bohrconv
#######################################################################################                        
        os.chdir(self.here)

        return (energy, gradients, scf_energy, nac, soc)

    def read_out_mol(self, pop="EPS"):
        """Read the output log file and return Mol"""
        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)
        out_mol = rf.mol_from_gauss(self.calc_name + ".log")
        os.chdir(self.here)

        return out_mol

    def read_hessian(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Get the Hessian matrix from a Gaussian output

        Returns
        ----------
        hessian : 3Natoms x 3Natoms array where Natoms is the amount of atoms in the
        QM region plus the atoms in the flexible QM' region. 
        """
        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)
        proc_fchk = subprocess.call("formchk -0 gck.chk gck.fchk", shell=True)
        hess_tmp = rf.read_hessian_g_fchk("gck.fchk")
        #truncate the hessian matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                hess = np.zeros((dim_flex,dim_flex))
            else:
                hess = np.zeros((len(positions),len(positions)))
        # Fix gradients units to Hartree/Angstrom
            hess[:len(positions),:len(positions)] = hess_tmp[:len(positions),:len(positions)]
        else:
            hess = hess_tmp[:len(positions),:len(positions)] 
        os.chdir(self.here)
        return hess

    def read_mu(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Read dipole moment and dipole derivatives from a G09/G16 output

        Returns
        ----------
        d_mu : array of dipole derivatives
        """

        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)
        proc_fchk = subprocess.call("formchk -0 gck.chk gck.fchk", shell=True)
        d_mu_tmp = rf.read_gauss_mu('gck.fchk')

        #truncate the dipole derivatives matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3 * natoms_flex)
                d_mu = np.zeros((dim_flex,3))
            else:
                d_mu = np.zeros((len(positions),3))
            d_mu[:len(positions),:3] = d_mu_tmp[:len(positions),:3]
        else:
            d_mu = d_mu_tmp[:len(positions),:3]

        os.chdir(self.here)
        return d_mu

    def read_nacs(self):
        """
        Read nonadiabatic coupling matrix elements from a TDDFT calculation from Gaussian16
        The matrix read is the not normalised.
        """
        gauss_path = os.path.join(self.here, self.calc_name)
        os.chdir(gauss_path)
        nacs = rf.read_gaussian_nacs("gck_ES.fchk")

        os.chdir(self.here)

        return nacs

class Gauss_CAS_calc(Calc):
    """
    Calculation with Gaussian 09 for CAS calculations
    """

    def run(self, atoms, nprocs=None):
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

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
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
    Calculation of TDDFT energy and gradients with Turbomole

    """
    def run(self, atoms, point_flex = None, nprocs = None, state=None, states=None, singlestate=0, nac_coupling=[], soc_coupling=[]):
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

        if state is not None and states is not None:
            ef.write_turbo_dyn("control", "control.temp", state, states, singlestate, nac_coupling, soc_coupling)

        # Run Turbomole
        proc = subprocess.Popen(
            "dscf > dscf.out && egrad > grad.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, point_flex = None, nprocs = None, state = None):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a normal modes 
        calculation.

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

        # Run normal modes calculation in the excited state
#        env = os.environ.copy()
#        env["state"] = state

        commands = ["actual -r",
            "dscf > dscf.out",
            "egrad > grad.out",
            "aoforce > force.out"
        ]

        for command in commands:
            result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True, env=env)

            if result.returncode != 0:print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            break

        os.chdir(self.here)

        return proc

    def read_out(self,
                 positions,
                 dyn_bool = False,
                 in_mol = None,
                 in_shell = None,
                 natoms_flex = None,
                 natoms = None,
                 state = None,
                 states = None,
                 mult = [],
                 singlestate = 0,
                 soc_coupling = []):
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

        nac = []
        soc = []

        if state is not None and states is not None:
            energy, gradients_b, scf_energy, nac, soc = rf.read_tb_dyn_tddft("grad.out",
                                                                              natoms,
                                                                              state,
                                                                              states,
                                                                              mult,
                                                                              singlestate,
                                                                              soc_coupling)

            # fix gradients units to Hartree/Angstrom
            gradients = gradients_b * bohrconv

        else:
            energy, gradients_b, scf_energy = rf.read_tb_grout("grad.out")
            # fix gradients units to Hartree/Angstrom
            gradients = gradients_b * bohrconv
            # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)
        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)
        return (energy, gradients, scf_energy, nac, soc)

    def read_mu(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Read dipole moment and dipole derivatives from a Turbomole output

        Returns
        ----------
        d_mu : array of dipole derivatives
        """

        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)
        d_mu_tmp = rf.read_turbo_mu('dipgrad')

        #truncate the dipole derivatives matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3 * natoms_flex)
                d_mu = np.zeros((dim_flex,3))
            else:
                d_mu = np.zeros((len(positions),3))
            d_mu[:len(positions),:3] = d_mu_tmp[:len(positions),:3]
        else:
            d_mu = d_mu_tmp[:len(positions),:3]

        os.chdir(self.here)

        return d_mu

class Turbo_calc_MP2(Calc):
    """
    Calculation with MP2 with Turbomole 7.0 and 7.6

    """
    def run(self, atoms, nprocs=None):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a SH dynamics,
        the placeholder &GRAD has to be included to specify for how 
        many states the gradient will be computed

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

        subprocess.call("cp control.temp control", shell=True)

#        turbo_redefine(atoms) FJH
        ef.write_coord(atoms)

        # Run Turbomole
        proc = subprocess.Popen(
            "jobex -level cc2 -c 1 > opt.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, point_flex = None, nprocs = None):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a normal modes 
        calculation.

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

        commands = ["actual -r",
            "dscf > dscf.out",
            "ricc2 > ricc2.out",
            "aoforce > force.out"
        ]

        for command in commands:
            result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True, env=env)

            if result.returncode != 0:print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            print(f"Error executing '{command}' in Turbo MP2: {result.stderr.decode()}")
            break

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyse a Turbomole job.last file while printing geometry updates

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

        energy, gradients_b, scf_energy = rf.read_tb_MP2_grout("job.last")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        subprocess.call("rm CC*", shell=True)

        os.chdir(self.here)
        return (energy, gradients, scf_energy)


class Turbo_calc(Calc):
    """
    Calculation with ADC2 or CC2 with Turbomole 7.0 and 7.6

    """

    def run(self, atoms, point_flex = None, nprocs = None, state=None, states=None, singlestate=0, nac_coupling=[], soc_coupling=[]):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a SH dynamics,
        the placeholder &GRAD has to be included to specify for how 
        many states the gradient will be computed

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

        if state is not None and states is not None:
            ef.write_turbo_dyn("control", "control.temp", state, states, singlestate, nac_coupling, soc_coupling)

        # Run Turbomole
        proc = subprocess.Popen(
            "dscf > dscf.out && ricc2 > ricc2.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, point_flex = None, nprocs = None, state = None):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a normal modes 
        calculation.

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

        # Run normal modes calculation in the excited state
#        env = os.environ.copy()
#        env["state"] = state

        commands = ["actual -r",
            "dscf > dscf.out",
            "ricc2 > ricc2.out",
            "aoforce > force.out"
        ]
     
#            "NumForce -ex $state -central -level cc2 > force.out"

        for command in commands:
            result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True, env=env)

            if result.returncode != 0:print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            break

        os.chdir(self.here)

        return proc

    def read_out(self,
                 positions,
                 dyn_bool = False,
                 in_mol = None,
                 in_shell = None,
                 natoms_flex = None,
                 natoms = None,
                 state = None,
                 states = None,
                 mult = [],
                 singlestate = 0,
                 soc_coupling = []):

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

        nac = []
        soc = []

        if state is not None and states is not None:
                        # read ouptut
            energy, gradients_b, scf_energy, nac, soc = rf.read_turbo_dyn("ricc2.out",
                                                                           natoms,
                                                                           state,
                                                                           states,
                                                                           mult,
                                                                           singlestate,
                                                                           soc_coupling)

            # fix gradients units to Hartree/Angstrom
            gradients = gradients_b * bohrconv

            # clean turbomole CC* files
#            subprocess.call("rm CC*", shell=True)
        else:
            energy, gradients_b, scf_energy = rf.read_ricc2("ricc2.out")
            # fix gradients units to Hartree/Angstrom
            gradients = gradients_b * bohrconv
            # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)
        # truncate gradients if too long
        gradients = gradients[:len(positions)]
 
        os.chdir(self.here)
        return (energy, gradients, scf_energy, nac, soc)


    def read_mu(self, positions, in_mol=None, in_shell=None, natoms_flex=None):        
        """
        Read dipole moment and dipole derivatives from a Turbomole output

        Returns
        ----------
        d_mu : array of dipole derivatives
        """

        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)
        d_mu_tmp = rf.read_turbo_mu('dipgrad')

        #truncate the dipole derivatives matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3 * natoms_flex)
                d_mu = np.zeros((dim_flex,3))
            else:
                d_mu = np.zeros((len(positions),3))
            d_mu[:len(positions),:3] = d_mu_tmp[:len(positions),:3]
        else:
            d_mu = d_mu_tmp[:len(positions),:3]

        os.chdir(self.here)

        return d_mu

    def save_checkpoint(self,
                        stp_iter = None,
                        Eini = None,
                        stp_ener = None,
                        stp_grad = None,
                        stp_grad_prev = None,
                        stp_vel = None,
                        stp_geom = None):
        """
        Save all the important info to restart a failed trajectory in a dynamics using OpenMolcas
 
        """
        turbo_path = os.path.join(self.here, self.calc_name)
        os.chdir(turbo_path)
        os.environ["curr_stp"] = str(stp_iter)
        subprocess.call("mkdir Step$curr_stp",shell=True)
        with open("dyn_restart","a") as check:
            check.write("%s" % "Checkpoint at step ")
            check.write("%s\n" % stp_iter)
            check.write("%s\n" % "----------")
            check.write("%s" % "Total Energy at step 1: ")
            check.write("%s\n" % Eini)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Electronic energies at current step")
            check.write("%s\n" % stp_ener)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Molecular Strucure")
            check.write("%s\n" % stp_geom)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Velocities")
            check.write("%s\n" % stp_vel)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Gradients")
            check.write("%s\n" % stp_grad)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Gradients at previous step")
            check.write("%s\n" % stp_grad_prev)
            check.write("%s\n" % "----------")
        check.closed
        subprocess.call("mv dyn_restart Step$curr_stp/",shell=True)
        subprocess.call("cp * Step$curr_stp/",shell=True)
        os.chdir(self.here)
        return

class Turbo_SCF_calc(Calc):
    """
    Calculation of SCF like DFT or HF with Turbomole

    """

    def run(self, atoms, nprocs=None):
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

        subprocess.call("cp control.temp control", shell=True)
        
        # Run Turbomole
        proc = subprocess.Popen(
            "dscf > dscf.out && grad > grad.out", stdout=FNULL, shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, point_flex = None, nprocs = None):
        """
        Write a Turbomole coord file and return a subprocess.Popen
        If the input file is going to be used for a normal modes 
        calculation.

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

        # Run normal modes calculation in the excited state

        commands = ["actual -r",
            "dscf > dscf.out",
            "grad > grad.out",
            "aoforce > force.out"
        ]

        for command in commands:
            result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, shell=True, env=env)

            if result.returncode != 0:print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            print(f"Error executing '{command}' in Turbo ADC2/CC2: {result.stderr.decode()}")
            break

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
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

        subprocess.call("mv gradient last_gradient", shell=True)
        subprocess.call("mv grad.out last_grad.out", shell=True)
        subprocess.call("mv dscf.out last_dscf.out", shell=True)
        os.chdir(self.here)
        return (energy, gradients, scf_energy)


class Molcas_calc(Calc):
    """
    Calculation with OpenMolcas
    """

    def __init__(self, calc_name_in=None, in_here=os.getcwd()):
        super().__init__(calc_name_in=calc_name_in, in_here=in_here)
        """
        It attempts to find molcas enviroment variables. Then, it wraps moclas calculation in mh/molcas
        molcas_path          molcas input files/template folder, in $PWD/mh
        molcas_project       molcas job name, default is molcas
        molcas_workdir       molcas scratch folder will be deleted after calculation done, default in /mh/molcas/molcas
        molcas_calcdir       molcas calculation folder in /mh/molcas, contains all output files
        """

        self.molcas_path = os.path.join(self.here, self.calc_name)
        self.molcas_calcdir = os.path.join(self.molcas_path, 'molcas')

        print("Molcas calculation run in %s" % (self.molcas_calcdir))

        try:
            self.molcas_project = os.environ["MOLCAS_PROJECT"]
            print("MOLCAS_PROJECT is found! Use job name: %s" % (self.molcas_project))
        except:
            self.molcas_project = 'molcas'
            print("MOLCAS_PROJECT is not set! Use default name: %s" % (self.molcas_project))
            os.environ["MOLCAS_PROJECT"] = self.molcas_project

        try:
            self.molcas_workdir = os.environ["MOLCAS_WORKDIR"]
            print("MOLCAS_WORKDIR is found! Write scratch in %s/%s" % (self.molcas_workdir, self.molcas_project))
        except:
            self.molcas_workdir = self.molcas_calcdir
            print("MOLCAS_WORKDIR is not set! Use default scratch: %s/%s" % (self.molcas_workdir, self.molcas_project))
            os.environ["MOLCAS_WORKDIR"] = '%s' % (self.molcas_workdir)

        self.molcas_scratch = os.path.join(self.molcas_workdir, self.molcas_project)

    def run(self, atoms, point_flex = None, nprocs=None, state=None, states=None, singlestate=0, nac_coupling=[], soc_coupling=[]): #FJH
        """
        Write a Molcas input file and return a subprocess.Popen

        Make sure the input file is called [name of calculation].input
        e.g. mh.input and the geometry file in Gateway is called geom.xyz

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Molcas
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method

        """
        molcas_path = os.path.join(self.here, self.calc_name)
        os.chdir(molcas_path)

        # Write a temporary geom file for molcas to read
        ef.write_xyz("geom.xyz", atoms)

        if point_flex is not None:
            if state is not None and states is not None:
                ef.write_molcas_free("molcas.input", self.calc_name + ".temp", state, states, singlestate, nac_coupling, soc_coupling, point_flex) 
            else: 
                ef.write_molcas("molcas.input", self.calc_name + ".temp", point_flex, freq = None) # The last bool turn the freq calc OFF
        else:
            if state is not None and states is not None:
                ef.write_molcas_free("molcas.input", self.calc_name + ".temp", state, states, singlestate, nac_coupling, soc_coupling, [])

        # Make molcas calculation folder and scratch if they do not exist 
        # as Molcas may not have the perssion to create it during calculation.
        if os.path.exists(self.molcas_calcdir) == False:
            subprocess.call("mkdir -p %s" % (self.molcas_calcdir), shell=True)

        if os.path.exists(self.molcas_scratch) == False:
            subprocess.call("mkdir -p %s" % (self.molcas_scratch), shell=True)

        # Copy input and coordiantes from $PWD/mh/ to $PWD/mh/molcas
        subprocess.call("cp -rf molcas.input %s/molcas.input" % (self.molcas_calcdir), shell=True)
        subprocess.call("cp -rf geom.xyz %s/geom.xyz" % (self.molcas_calcdir), shell=True)

        # Copy orbital from $PWD/mh to $PWD/mh/molcas
        # Use StrOrb if RasOrb does not exist
        # return error if StrOrb and RasOrb do not exist
        if os.path.exists('%s.StrOrb' % (self.molcas_project)) == True and os.path.exists('%s.RasOrb' % (self.molcas_project)) == False:
            subprocess.call("cp -rf %s.StrOrb %s/%s.RasOrb" % (self.molcas_project, self.molcas_calcdir, self.molcas_project), shell=True)
            print('Molcas found orbital file %s.StrOrb' % (self.molcas_project))

        elif os.path.exists('%s.StrOrb' % (self.molcas_project)) == False and os.path.exists('%s.RasOrb' % (self.molcas_project)) == False:
            import sys
            sys.exit('FileNotFoundError\n  Molcas is looking for guess orbital %s.StrOrb or %s.RasOrb' % (self.molcas_project, self.molcas_project))

        elif os.path.exists('%s.RasOrb' % (self.molcas_project)) == True:
            print('Molcas found orbital file %s.RasOrb' % (self.molcas_project))
            subprocess.call("cp -rf %s.RasOrb %s/%s.RasOrb" % (self.molcas_project, self.molcas_calcdir, self.molcas_project), shell=True)

        # run molcas calculation in molcas_calcdir
        # add '-b 1' to print output on-the-fly
        os.chdir(self.molcas_calcdir)
        os.environ["np"] = nprocs

        if state is not None and states is not None:
            proc = subprocess.Popen(
            "pymolcas -nt $np molcas.input -f -b 1", shell=True)
        else:
#            os.environ["np"] = nprocs
            proc = subprocess.Popen(
                "pymolcas -nt $np molcas.input -f -b 1", shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self,atoms, point_flex = None, nprocs=None):
        """
        Run a normal modes calculation in Molcas
        """
        molcas_path = os.path.join(self.here, self.calc_name)
        os.chdir(molcas_path)

        # Write a temporary geom file for molcas to read
        ef.write_xyz("geom.xyz", atoms)
  
        freq = True
        if point_flex is not None:
            ef.write_molcas("molcas.input", self.calc_name + ".temp", point_flex, freq) # The last value turn the freq calc ON
        else:
            ef.write_molcas("molcas.input", self.calc_name + ".temp", [], freq)

        # Make molcas calculation folder and scratch if they do not exist 
        # as Molcas may not have the perssion to create it during calculation.
        if os.path.exists(self.molcas_calcdir) == False:
            subprocess.call("mkdir -p %s" % (self.molcas_calcdir), shell=True)

        if os.path.exists(self.molcas_scratch) == False:
            subprocess.call("mkdir -p %s" % (self.molcas_scratch), shell=True)

        # Copy input and coordiantes from $PWD/mh/ to $PWD/mh/molcas
        subprocess.call("cp -rf molcas.input %s/molcas.input" % (self.molcas_calcdir), shell=True)
        subprocess.call("cp -rf geom.xyz %s/geom.xyz" % (self.molcas_calcdir), shell=True)

        # Copy orbital from $PWD/mh to $PWD/mh/molcas
        # Use StrOrb if RasOrb does not exist
        # return error if StrOrb and RasOrb do not exist
        if os.path.exists('%s.StrOrb' % (self.molcas_project)) == True and os.path.exists('%s.RasOrb' % (self.molcas_project)) == False:
            subprocess.call("cp -rf %s.StrOrb %s/%s.RasOrb" % (self.molcas_project, self.molcas_calcdir, self.molcas_project), shell=True)
            print('Molcas found orbital file %s.StrOrb' % (self.molcas_project))

        elif os.path.exists('%s.StrOrb' % (self.molcas_project)) == False and os.path.exists('%s.RasOrb' % (self.molcas_project)) == False:
            import sys
            sys.exit('FileNotFoundError\n  Molcas is looking for guess orbital %s.StrOrb or %s.RasOrb' % (self.molcas_project, self.molcas_project))

        elif os.path.exists('%s.RasOrb' % (self.molcas_project)) == True:
            print('Molcas found orbital file %s.RasOrb' % (self.molcas_project))
            subprocess.call("cp -rf %s.RasOrb %s/%s.RasOrb" % (self.molcas_project, self.molcas_calcdir, self.molcas_project), shell=True)

        # run molcas calculation in molcas_calcdir
        # add '-b 1' to print output on-the-fly
        os.chdir(self.molcas_calcdir)

        #if state is not None and states is not None:
        proc = subprocess.Popen(
            "pymolcas -nt molcas.input -f -b 1", shell=True)

        #os.environ["np"] = nprocs
        #proc = subprocess.Popen(
        #    "pymolcas -np $np molcas.input -f -b 1", shell=True)

        return proc

    def read_out(self,
                 positions,
                 dyn_bool = False,
                 in_mol = None,
                 in_shell = None,
                 natoms_flex = None,
                 natoms = None,
                 state = None,
                 states = None,
                 mult = [],
                 singlestate = 0,
                 soc_coupling = []):
        """
        Analyse a Molcas .log file while printing geometry updates

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

        os.chdir(self.molcas_path)

        #initialize nac and soc, which will be used when dyn_bool == True
        nac = []
        soc = []

        subprocess.call("cp -rf %s/%s.RasOrb %s.RasOrb" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)
        subprocess.call("cp -rf %s/molcas.log molcas.log" % (self.molcas_calcdir), shell=True)
        subprocess.call("cp -rf %s/%s.rasscf.h5 %s.rasscf.h5" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)
        subprocess.call("cp -rf %s/%s.rasscf.molden %s.rasscf.molden" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)

        if state is not None and states is not None:

            # read ouptut
            energy, gradients_b, scf_energy, nac, soc = rf.read_molcas_ext("molcas.log",
                                                                           natoms,
                                                                           state,
                                                                           states,
                                                                           mult,
                                                                           singlestate,
                                                                           soc_coupling)

            # fix gradients units to Hartree/Angstrom
            gradients = gradients_b * bohrconv


            # clean molcas scratch
            # note that the removal of molcas_calcdir will be controlled by fro_dyn 
            subprocess.call("rm -rf %s" % (self.molcas_scratch), shell=True)

        else:
            energy, gradients_b, scf_energy = rf.read_molcas("molcas.log")
            # fix gradients units to Hartree/Angstrom
#            gradients = gradients_b * bohrconv
            # update the geometry log
            if in_mol != None:
                self.update_geom(positions, in_mol, in_shell)

            # truncate gradients if too long and fix gradients units to Hartree/Angstrom
            if natoms_flex is not None:
                if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                    dim_flex = int(len(positions) + 3. * natoms_flex)
                    gradients = np.zeros(dim_flex)
                else:
                    gradients = np.zeros(len(positions))
                gradients[:len(positions)] = gradients_b[:len(positions)] * bohrconv
            else:
                gradients = gradients_b[:len(positions)] * bohrconv

 #           if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
 #               dim_flex = int(len(positions) + 3. * natoms_flex)
 #               gradients = np.zeros(dim_flex)
 #           else:
 #               gradients = np.zeros(len(positions))
            # Fix gradients units to Hartree/Angstrom
            # truncate gradients if too long
 #           gradients[:len(positions)] = gradients_b[:len(positions)] * bohrconv


#            dim_flex = int(len(positions) + 3. * natoms_flex)
#            gradients = np.zeros(dim_flex)
            # truncate gradients if too long
#            gradients[:len(positions)] = gradients_b[:len(positions)] * bohrconv
#            gradients = gradients[:len(positions)]
        os.chdir(self.here)

        return (energy, gradients, scf_energy, nac, soc)

    def read_hessian(self):
        """
        Get the Hessian matrix from a Molcas output
        In Molcas, the hessian is stored in the 
        "name".slapaf.h5 file

        Returns
        ----------
        hessian : 3Natoms x 3Natoms array where Natoms is the amount of atoms in the
        QM region plus the atoms in the flexible QM' region.
        """

        os.chdir(self.molcas_path)

        subprocess.call("cp -rf %s/%s.RasOrb %s.RasOrb" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)
        subprocess.call("cp -rf %s/molcas.log molcas.log" % (self.molcas_calcdir), shell=True)
        subprocess.call("cp -rf %s/%s.slapaf.h5 %s.slapaf.h5" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)
        subprocess.call("cp -rf %s/%s.rasscf.h5 %s.rasscf.h5" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)
        subprocess.call("cp -rf %s/%s.rasscf.molden %s.rasscf.molden" % (self.molcas_calcdir, self.molcas_project, self.molcas_project), shell=True)

        hess_tmp = rf.read_hessian_molcas("molcas.slapaf.h5")
        #truncate the hessian matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                hess = np.zeros((dim_flex,dim_flex))
            else:
                hess = np.zeros((len(positions),len(positions)))
            # Fix gradients units to Hartree/Angstrom
            hess[:len(positions),:len(positions)] = hess_tmp[:len(positions),:len(positions)]
        else:
            hess = hess_tmp[:len(positions),:len(positions)]
        os.chdir(self.here)
        return hess

    def read_nacs(self):
        """
        Read nonadiabatic coupling matrix elements from a CASSCF calculation
        The matrix read is the not normalised one. The option inclusing CSF 
        is not implemented yet
        """
        molcas_path = os.path.join(self.here, self.calc_name)
        os.chdir(molcas_path)
        nacs = rf.read_molcas_nacs("molcas.log")

        os.chdir(self.here)

        return nacs
 
    def save_checkpoint(self,
                        stp_iter = None,
                        Eini = None,
                        stp_ener = None,
                        stp_grad = None,
                        stp_grad_prev = None,
                        stp_vel = None,
                        stp_geom = None):
        """
        Save all the important info to restart a failed trajectory in a dynamics using OpenMolcas
 
        """
        molcas_path = os.path.join(self.here, self.calc_name)
        os.chdir(molcas_path)
        os.environ["curr_stp"] = str(stp_iter)
        subprocess.call("mkdir Step$curr_stp",shell=True)
        subprocess.call("cp molcas.RasOrb* Step$curr_stp/",shell=True)
        subprocess.call("cp molcas.rasscf.* Step$curr_stp/",shell=True)
        with open("dyn_restart","a") as check:
            check.write("%s" % "Checkpoint at step ")
            check.write("%s\n" % stp_iter)
            check.write("%s\n" % "----------")
            check.write("%s" % "Total Energy at step 1: ")
            check.write("%s\n" % Eini)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Electronic energies at current step")
            check.write("%s\n" % stp_ener)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Molecular Strucure")
            check.write("%s\n" % stp_geom)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Velocities")
            check.write("%s\n" % stp_vel)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Gradients")
            check.write("%s\n" % stp_grad)
            check.write("%s\n" % "----------")
            check.write("%s\n" % "Gradients at previous step")
            check.write("%s\n" % stp_grad_prev)
            check.write("%s\n" % "----------")
        check.closed
        subprocess.call("mv dyn_restart Step$curr_stp/",shell=True)
        os.chdir(self.here)
        return

    def try_reord_act_spc(self,alter_orbs):
        moclas_path = os.path.join(self.here,self.calc_name)
        os.chdir(molcas_path)

        file_name = self.calc_name + ".temp"
        with open(file_name) as temp_file:
            temp_content = temp_file.readlines()
        out_file = open("mh.temp","w")
        for line in temp_content:
            if "XXX_ALTER_ORBS_XXX" in line:
                out_file.write(line.replace("XXX_ALTER_ORBS_XXX",str(alter_orbs)))
            else:
                out_file.write(line)
        out_file.close()
        os.chdir(self.here)
        return
          
class xtb_calc(Calc):
    """
    Calculation of energy and gradients with GFN2-xTB

    """

    def run(self, atoms, point_flex = None, nprocs=None):
        """
        Write a GFN2-xTB coordinate file and return a subprocess.Popen
        
        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with GFN2-xTB
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)

        # Write temporary geometry file in XYZ format
        if point_flex is not None:
            ef.write_xtb("geom.xyz", atoms,
                           point_flex, self.calc_name + ".temp")
        else:
            ef.write_xtb("geom.xyz", atoms,
                           [], self.calc_name + ".temp")

        if os.path.isfile("xtb.input"):
            xtb_run_string = "xtb -I xtb.input geom.xyz --grad --iterations 2000 --norestart > xtb.out"
        else:
            xtb_run_string = "xtb geom.xyz --grad --iterations 2000 --acc 10 --norestart > xtb.out"

        proc = subprocess.Popen(
            xtb_run_string, shell=True)

        os.chdir(self.here)
        
        return proc
 
    def run_freq(self,atoms, point_flex = None, nprocs=None):
        """
        Write a GFN2-xTB coordinate file and return a subprocess.Popen
        to compute the Hessian
        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with GFN2-xTB
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)

        # Write temporary geometry file in XYZ format
        if point_flex is not None:
            ef.write_xtb("geom.xyz", atoms,
                           point_flex, self.calc_name + ".temp")
        else:
            ef.write_xtb("geom.xyz", atoms,
                           [], self.calc_name + ".temp")

        if os.path.isfile("xtb.input"):
            xtb_run_string = "xtb -I xtb.input geom.xyz --etemp 100 --iterations 2000 --hess --norestart > xtb.out"
        else:
            xtb_run_string = "xtb geom.xyz --etemp 100 --iterations 2000 --hess --norestart > xtb.out"

        proc = subprocess.Popen(
            xtb_run_string, shell=True)

        os.chdir(self.here)

        return proc   
        

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyze a GFN2-xTB gradients file  while printing
        geometry updates

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
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)
        
        energy, gradients_bohr = rf.read_xtb("gradient")
        scf_energy = energy
        # Fix gradients units to Hartree/Angstrom
#        gradients = gradients_bohr * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)


        # truncate gradients if too long and fix gradients units to Hartree/Angstrom
        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                gradients = np.zeros(dim_flex)
            else:
                gradients = np.zeros(len(positions))
            gradients[:len(positions)] = gradients_bohr[:len(positions)] * bohrconv
        else:
            gradients = gradients_bohr[:len(positions)] * bohrconv
        # Fix gradients units to Hartree/Angstrom
#        gradients[:len(positions)] = gradients_bohr[:len(positions)] * bohrconv

        os.chdir(self.here)
 
        return (energy, gradients, scf_energy)

    def read_charges(self):
        """
        Get the atomic charges of the whole system
        
        Returns
        ----------
        charges : array of atom charges
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)
        charges = rf.read_xtb_charges("charges")
        os.chdir(self.here)
        return charges
   
    def read_hessian(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Get the Hessian matrix from a xTB output

        Returns
        ----------
        hess : 3Natoms x 3Natoms array where Natoms is the amount of atoms in the
        QM region plus the atoms in the flexible QM' region.
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)
        hess_tmp = rf.read_hessian_xtb("hessian")
        #truncate the hessian matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                hess = np.zeros((dim_flex,dim_flex))
            else:
                hessian = np.zeros((len(positions),len(positions)))
            # Fix gradients units to Hartree/Angstrom
            hess[:len(positions),:len(positions)] = hess_tmp[:len(positions),:len(positions)]
        else:
            hess = hess_tmp[:len(positions),:len(positions)]

        os.chdir(self.here)
        return hess

class xtb_calc_gfnff(Calc):
    """
    Calculation of energy and gradients with GFN2-FF from xTB

    """

    def run(self, atoms, point_flex = None, nprocs=None):
        """
        Write a xTB GFN-FF coordinate file and return a subprocess.Popen
        
        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with xTB GFN-FF
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)

        # Write temporary geometry file in XYZ format
        if point_flex is not None:
            ef.write_xtb("geom.xyz", atoms,
                           point_flex, self.calc_name + ".temp")
        else:
            ef.write_xtb("geom.xyz", atoms,
                           [], self.calc_name + ".temp")

        if os.path.isfile("xtb.input"):
            xtb_run_string = "xtb -I xtb.input geom.xyz --grad --gfnff --iterations 2000 --norestart > xtb.out"
        else:
            xtb_run_string = "xtb geom.xyz --grad --gfnff --iterations 2000 --acc 10 --norestart > xtb.out"

        proc = subprocess.Popen(
            xtb_run_string, shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self,atoms, point_flex = None, nprocs=None):
        """
        Write a xTB GFN-FF coordinate file and return a subprocess.Popen
        to compute the Hessian
        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with xTB GFN-FF
        Returns
        -------
        proc : subprocess.Popen object
            the object should have a .wait() method
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)

        # Write temporary geometry file in XYZ format
        if point_flex is not None:
            ef.write_xtb("geom.xyz", atoms,
                           point_flex, self.calc_name + ".temp")
        else:
            ef.write_xtb("geom.xyz", atoms,
                           [], self.calc_name + ".temp")

        if os.path.isfile("xtb.input"):
            xtb_run_string = "xtb -I xtb.input geom.xyz --gfnff --etemp 100 --iterations 2000 --hess --norestart > xtb.out"
        else:
            xtb_run_string = "xtb geom.xyz --gfnff --etemp 100 --iterations 2000 --hess --norestart > xtb.out"

        proc = subprocess.Popen(
            xtb_run_string, shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyze a xTB GFN-FF gradients file  while printing
        geometry updates

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
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)

        energy, gradients_bohr = rf.read_xtb("gradient")
        scf_energy = energy
        # Fix gradients units to Hartree/Angstrom
#        gradients = gradients_bohr * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long and fix gradients units to Hartree/Angstrom
        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                gradients = np.zeros(dim_flex)
            else:
                gradients = np.zeros(len(positions))
            gradients[:len(positions)] = gradients_bohr[:len(positions)] * bohrconv
        else:
            gradients = gradients_bohr[:len(positions)] * bohrconv
        # Fix gradients units to Hartree/Angstrom
#        gradients[:len(positions)] = gradients_bohr[:len(positions)] * bohrconv

        os.chdir(self.here)

        return (energy, gradients, scf_energy)

    def read_charges(self):
        """
        Get the atomic charges of the whole system
        
        Returns
        ----------
        charges : array of atom charges
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)
        charges = rf.read_xtb_charges("gfnff_charges")
        os.chdir(self.here)
        return charges

    def read_hessian(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Get the Hessian matrix from a xTB GFN-FF output

        Returns
        ----------
        hess : 3Natoms x 3Natoms array where Natoms is the amount of atoms in the
        QM region plus the atoms in the flexible QM' region.
        """
        xtb_path = os.path.join(self.here, self.calc_name)
        os.chdir(xtb_path)
        hess_tmp = rf.read_hessian_xtb("hessian")
        #truncate the hessian matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                hess = np.zeros((dim_flex,dim_flex))
            else:
                hess = np.zeros((len(positions),len(positions)))
            # Fix gradients units to Hartree/Angstrom
            hess[:len(positions),:len(positions)] = hess_tmp[:len(positions),:len(positions)]
        else:
            hess = hess_tmp[:len(positions),:len(positions)]

        os.chdir(self.here)
        return hess


class Qchem(Calc):
    """
    Constrained DFT calculation performed with QChem
    """
    def run(self, atoms,nprocs):
        """
        Write a QChem input file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with QChem
        Returns
        -------
        proc : subprocess.Popen object
        """

        qchem_path = os.path.join(self.here, self.calc_name)
        os.chdir(qchem_path)
        # Writes modified qchem input
        ef.write_qchem(self.calc_name + ".in", atoms,
                       "qchem.temp")
        os.environ["np"] = nprocs
        proc = subprocess.Popen(
            "qchem -nt $np " + self.calc_name + ".in" + " " + self.calc_name + ".out", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyse a QChem .out file while printing geometry updates

        This is a modified function method from the Gaussian Calc class

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
            Energy calculated by QChem in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        qchem_path = os.path.join(self.here, self.calc_name)
        os.chdir(qchem_path)

        # energies are in Hartree
        # gradients are in Hartree/Bohr
        energy, gradients_b, scf_energy = rf.read_qchem_out(self.calc_name + ".out")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)

        return (energy, gradients, scf_energy)

class nwchem_calc_DFT(Calc):
    """
    DFT, CDFT and TDDFT calculations performed with NWChem
    """

    def run(self, atoms, nprocs):
        """
        Write a NWChem input file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with QChem
        Returns
        -------
        proc : subprocess.Popen object
        """

        nwchem_path = os.path.join(self.here, self.calc_name)
        os.chdir(nwchem_path)
        # Writes modified qchem input
        ef.write_nwchem(self.calc_name + ".nw", atoms,
                       "nw.temp")
        os.environ["np"] = nprocs
        proc = subprocess.Popen(
            "mpirun -np $np nwchem " + self.calc_name + ".nw" + " " + ">" + self.calc_name + ".out", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyse a NWChem .out file while printing geometry updates

        This is a modified function method from the Gaussian Calc class

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
            Energy calculated by NWChem in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        nwchem_path = os.path.join(self.here, self.calc_name)
        os.chdir(nwchem_path)

        # energies are in Hartree
        # gradients are in Hartree/Bohr
        energy, gradients_b, scf_energy = rf.read_nwchem_DFT_out(self.calc_name + ".out")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.chdir(self.here)

        return (energy, gradients, scf_energy)

class fomo_ci_calc(Calc):
    """
    MNDO, AM1, PM3, PM6 and FOMO-CI with all the previous Hamiltonians
    performed with MOPAC software
    """

    def run(self, atoms, nprocs, at_reparam=None):
        """
        Write a MOPAC input file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with QChem
        Returns
        -------
        proc : subprocess.Popen object
        """

        mopac_path = os.path.join(self.here, self.calc_name)
        os.chdir(mopac_path)
 

        mol = Mol(atoms)

        if self.calc_name == 'rl':
            # Write modified mopac inputs
            ef.write_mopac(self.calc_name + ".dat", atoms,"mopac.temp")
        else:
        # Writes modified mopac inputs
            if at_reparam is not None:
               for k in at_reparam:
                   k -= 1
                   atoms[k].elem = atoms[k].elem+"w" 
            ef.write_mopac(self.calc_name + ".dat", atoms,"mopac.temp")
            ef.write_tinker_xyz(self.calc_name + "_tnk.xyz", atoms, "mopac_tnk.temp")

        os.environ["np"] = nprocs
        proc = subprocess.Popen(
            "mpirun -np $np mopac2002.x " + self.calc_name + ".dat" , shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyse a MOPAC-FOMO-CI .out file while printing geometry updates

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
            Energy calculated by MOPAC in a.u.
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in a.u.
            NB: if multiplicity is set to triplets, the scf_energy is the
                value for T1 instead of S0.
        """
        mopac_path = os.path.join(self.here, self.calc_name)
        os.chdir(mopac_path)

        # energies are in Hartree
        # gradients are in Hartree/Bohr
        energy, gradients_b, scf_energy = rf.mopac_fomo_ci_out(self.calc_name + ".dat.out")
        # fix gradients units to Hartree/Angstrom
        gradients = gradients_b * bohrconv
        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        # truncate gradients if too long
        gradients = gradients[:len(positions)]

        os.environ["calc_name"] = self.calc_name
        subprocess.call("rm $calc_name.arc*", shell=True)

        os.chdir(self.here)

        return (energy, gradients, scf_energy)

    def read_nacs(self):
        """
        Yet to be implemented

        """
        pass

class Orca_calc(Calc):
    """
    DFT, TDDFT, SF-DFT and CASSCF calculations computed with Orca
    """
    def run(self, atoms, points_flex = None, nprocs=None):
        """
        Write a Orca input file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Orca
        Returns
        -------
        proc : subprocess.Popen object
        """

        orca_path = os.path.join(self.here, self.calc_name)
        os.chdir(orca_path)
        ef.write_orca(self.calc_name + ".inp", atoms,
                       "mh.temp")
        if points_flex is not None:
            ef.write_orca_charges(("charges.pc", point_flex))
        os.environ["np"] = nprocs
        proc = subprocess.Popen(
            "orca " + self.calc_name + ".inp" + " > " + self.calc_name + ".out", shell=True)

        os.chdir(self.here)

        return proc

    def run_freq(self, atoms, points_flex = None, nprocs=None):
        """
        Write a Orca input file and return a subprocess.Popen

        Parameters
        ----------
        atoms : list of Atom objects
            Atoms to be calculated with Orca
        Returns
        -------
        proc : subprocess.Popen object
        """

        orca_path = os.path.join(self.here, self.calc_name)
        os.chdir(orca_path)
        ef.write_orca(self.calc_name + ".inp", atoms,
                       "mh.temp")
        if points_flex is not None:
            ef.write_orca_charges(("charges.pc", point_flex))
        os.environ["np"] = nprocs
        proc = subprocess.Popen(
            "orca " + self.calc_name + ".inp" + " > " + self.calc_name + ".out", shell=True)

        os.chdir(self.here)

        return proc

    def read_out(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Analyse a Orca.out file while printing geometry updates

        This is a modified function method from the Gaussian Calc class

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
            Energy calculated by Orca in Hartree
        gradients : list of floats
            The gradients in form x1,y1,z1,x2,y2,z2 etc. in Hartree/Angstrom
        scf_energy : float
            The ground state energy in Hartree

        """
        orca_path = os.path.join(self.here, self.calc_name)
        os.chdir(orca_path)

        # energies are in Hartree
        # gradients are in Hartree/Bohr
        energy, gradients_bohr, scf_energy = rf.read_orca_out(self.calc_name + ".out")

        # truncate gradients if too long and fix gradients units to Hartree/Angstrom
        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3 * natoms_flex)
                gradients = np.zeros(dim_flex)
            else:
                gradients = np.zeros(len(positions))
            # Fix gradients units to Hartree/Angstrom
            gradients[:len(positions)] = gradients_bohr[:len(positions)] * bohrconv
        else:
            # Fix gradients units to Hartree/Angstrom
            gradients = gradients_bohr[:len(positions)] * bohrconv

        # update the geometry log
        if in_mol != None:
            self.update_geom(positions, in_mol, in_shell)

        os.chdir(self.here)

        return (energy, gradients, scf_energy)

    def read_hessian(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Get the Hessian matrix from an Orca output

        Returns
        ----------
        hessian : 3Natoms x 3Natoms array where Natoms is the amount of atoms in the
        QM region plus the atoms in the flexible QM' region. 
        """
        orca_path = os.path.join(self.here, self.calc_name)
        os.chdir(orca_path)
        hess_tmp = rf.read_hessian_orca(self.calc_name + ".hess")
        #truncate the hessian matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3. * natoms_flex)
                hess = np.zeros((dim_flex,dim_flex))
            else:
                hess = np.zeros((len(positions),len(positions)))
        # Fix gradients units to Hartree/Angstrom
            hess[:len(positions),:len(positions)] = hess_tmp[:len(positions),:len(positions)]
        else:
            hess = hess_tmp[:len(positions),:len(positions)]
        os.chdir(self.here)
        return hess

    def read_mu(self, positions, in_mol=None, in_shell=None, natoms_flex=None):
        """
        Read dipole moment and dipole derivatives from a Turbomole output
       
        Returns
        ----------
        d_mu : array of dipole derivatives
        """

        orca_path = os.path.join(self.here, self.calc_name)
        os.chdir(orca_path)
        d_mu_tmp = rf.read_orca_mu(self.calc_name + '.hess')

        #truncate the dipole derivatives matrix if it is too long

        if natoms_flex is not None:
            if int(len(positions)) < int(3*natoms_flex): # CHANGE THIS AWFULNESS PLEASE!
                dim_flex = int(len(positions) + 3 * natoms_flex)
                d_mu = np.zeros((dim_flex,3))
            else:
                d_mu = np.zeros((len(positions),3))
            d_mu[:len(positions),:3] = d_mu_tmp[:len(positions),:3]
        else:
            d_mu = d_mu_tmp[:len(positions),:3]

        os.chdir(self.here)
        return
