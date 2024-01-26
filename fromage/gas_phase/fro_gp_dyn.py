#!/usr/bin/env python
## Extension library for fromage for ab initio Molecular Dynamics
## based largely (with permission) on the PyRAI2MD package.
## 
## Developers
## Dr. Jordan Cox
## Dr. Jingbai Li
## Dr. Federico Hernandez

import time,datetime,os
import numpy as np
from fromage.dynamics.periodic_table import Element
from fromage.dynamics.verlet import NoseHoover, VerletI, VerletII
from fromage.dynamics.surfacehopping import FSSH, GSH, NOSH
from fromage.dynamics.tools import Printcoord, NACpairs
from fromage.utils import calc
from fromage.utils import array_operations as ao
from fromage.io import read_file as rf
from fromage.utils.atom import Atom
from fromage.io.parse_config_file import bool_cast

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


###################################################
############# Input Parameters Parser #############
###################################################

def initTrajParams(symbols, init_pos, init_vel, settings, curr_step=None, Eini=None, prev_grad=None):
    """
    Parse input parameters and return dict of settings

    Parameters
    ----------
    symbols : List<str>
        Atomic symbols for all atoms in high layer
    init_pos : List<float>
        Initial positions of atoms in high layer as [X1, Y1, Z1, X2, Y2, Z2...]
        units of coordinates are Angstroms
    init_vel : List<float>
        Initial XYZ components of atomic velocities as [[Vx1,Vy1,Vz1],[Vx2,Vy2,Vz2]...]
        units of velocity are Bohr/a.u. of time
    settings : Dict
        Unfiltered dict read from fromage.in file

    Returns
    -------
    ip : Dict of trajectory settings

    """
    ip = {}
    ip["types"] = symbols
    ip["R"] = np.array(init_pos)
    ip["M"] = np.array([getMass(x) for x in symbols])
    ip["V"] = np.array(init_vel)
    ip["State"] = settings["init_state"]
    ip["natoms"] = len(ip.get("types"))
    ip["stepsize"] = settings["step"]
    ip["T"] = settings["time"]
    ip["HopType"] = settings["hop_method"]
    ip["low_level"] = settings["low_level"]
    ip["high_level"] = settings["high_level"]
    ip["out_file"] = settings["out_file"]
    if "temp" in settings.keys():
        ip["temp"] = settings["temp"]
    if "check" in settings.keys():
        ip["check"] = settings["check"]
    if "verbose" in settings.keys():
        ip["verbose"] = settings["verbose"]
    if "hop_thresh" in settings.keys():
        ip["gap"] = settings["hop_thresh"]
    if "isc_tresh" in settings.keys():
        ip["gapsoc"] = settings["isc_tresh"]
    if "adj_mom" in settings.keys():
        ip["adjust"] = settings["adj_mom"]
    if "maxhop" in settings.keys():
        ip["maxhop"] = settings["maxhop"]
    if "decoherence" in settings.keys():
        ip["decoherence"] = settings["decoherence"]
    if "substeps" in settings.keys():
        ip["substeps"] = settings["substeps"]
    if "integrate" in settings.keys():
        ip["integrate"] = settings["integrate"]
    if "reflect" in settings.keys():
        ip["reflect"] = settings["reflect"]
    if "chk_stp" in settings.keys():
        ip["chk_stp"] = settings["chk_stp"]
    if "alter_orbs" in settings.keys():
        ip["alter_orbs"] = settings["alter_orbs"]
    if "e_cons" in settings.keys():
        ip["e_cons"] = settings["e_cons"]
    if "dyn_restart" in settings.keys():
        ip["dyn_restart"] = bool_cast(settings["dyn_restart"])
    if "nprocs" in settings.keys():
        ip["nprocs"] = settings["nprocs"]
    if "natoms_flex" in settings.keys():
        ip["natoms_flex"] = settings["natoms_flex"]
    if prev_grad is not None:
        ip["Gp"] = np.array(prev_grad)
    if "curr_step" != None:
        ip["curr_step"] = curr_step
    if "Eini" != None:
        ip["Eini"] = Eini
    if "nactype" in settings.keys():
        ip["nactype"] = str(settings["nactype"])
    if "singlestate" in settings.keys():
        ip["singlestate"] = int(settings["singlestate"])
    else:
        ip["singlestate"] = 0
    if "stop_traj" in settings.keys():
        ip["stop_traj"] = float(settings["stop_traj"])
    else:
        ip["stop_traj"] = None
    ## prepare spin states and couplings
    ## spin and states are list
    spin = [0]
    states = [2]
    statemult = []

    if "spin" in settings.keys():
        spin = [int(x) for x in settings["spin"]]
    if "states" in settings.keys():
        states = [int(x) for x in settings["states"]]

    mult = []
    nstates = int(np.sum(states))

    for n, s in enumerate(states):
        ms = int(spin[n] * 2 + 1)
        mult.append(ms)
        for m in range(s):
            statemult.append(ms)

    coupling = []
    nac_coupling = []
    soc_coupling = []

    if "coupling" in settings.keys():
        coupling = [int(x) for x in settings["coupling"]]
        coupling = np.array(coupling).reshape((-1, 2))

    for pair in coupling:
        s1, s2 = pair
        s1 -= 1
        s2 -= 1
        if statemult[s1] != statemult[s2]:
            soc_coupling.append(sorted([s1, s2]))
        else:
            nac_coupling.append(sorted([s1, s2]))

    ip["states"] = states
    ip["nstates"] = nstates
    ip["mult"] = mult
    ip["statemult"] = statemult
    ip["nac_coupling"] = nac_coupling
    ip["soc_coupling"] = soc_coupling

    return ip

####################################################
########### Energy Function for Dynamics ###########
####################################################

def dynamics_sequence(traj):
    """
    Run QM calculations in parallel and write and return results

    This function is a modified version of the sequence() function in
    fro_run.py which computes multistate energies and gradients, and 
    nonadiabatic couplings, for use with nonadiabatic dynamics simulations.

    Parameters
    ----------
    in_pos : list<float>
        Input coordinates in 1D array form
    mol_atoms : list<Atom>
        List of Atom objects that acts as a template for high layer
    state : int
        Index of current electronic state
    nstates : int
        Number of electronic states in state averaging scheme
    high_level : str
        Program to use for high level calculation
    read_nacs : bool
        Flag for reading NACs from mh calculation, only for nactype == 'nac'

    Returns
    -------
    en_out : list<float>
        List of electronic energies for each state in Hartree
    gr_out : list<float>
        List of gradients of en_out in Hartree/Angstrom
    nac_out : list<float>
        List of nonadiabatic couplings between electronic states
    soc_out : list<float>
        List of spin-orbit coupling between spin states
    """
    SH_methods = ['molcas','turbomole','turbomole_tddft','qchem','gaussian']
    # Read parameters from Trajectory object
    in_pos = atoms_to_fromage(traj.R.copy())
    mol_atoms = traj.mol_atoms
    state = traj.state
    states = traj.states
    mult = traj.mult
    statemult = traj.statemult
    nstates = traj.nstates
    high_level = traj.high
    out_file = open(traj.out_file,'a+')
    natoms = traj.natoms
    singlestate = traj.singlestate
    nactype = traj.nactype
    nac_coupling = traj.nac_coupling
    soc_coupling = traj.soc_coupling
    nprocs = traj.nprocs
    stop_traj = traj.stop_traj
    natoms_flex = None
    
    # Check the method selected for the high_level is supported for SH-dynamics
    if high_level in SH_methods:
        pass
    else:
        out_file.write(" The method %s is not implemented for SH dynamics\n" % (high_level))
        out_file.write("The job is dying now :-( ")
        import sys
        sys.exit()
    # initialise calculation objects
    mh = calc.setup_calc("mh", high_level)

    # Run the calculations as subprocesses in parallel
    calcs = []

    pass_nac = []
    if nactype == 'nac':
        pass_nac = nac_coupling

    mh_proc = mh.run(ao.array2atom(mol_atoms, in_pos),
                     nprocs = nprocs,
                     state = state,
                     states = states,
                     singlestate = singlestate,
                     nac_coupling = pass_nac,
                     soc_coupling = soc_coupling)
    calcs.append(mh_proc)

    ## Wait until all parallel calculations are finished
    for proc in calcs:
        proc.communicate()

     # read results. Each x_en_gr is a tuple (energy,gradients,scf_energy)

    if high_level == "gaussian_cas":
        mh_en_gr = mh.read_out(in_pos, natoms_r2 = natoms_r2)[0:3]
    elif high_level in SH_methods:
#    elif high_level == "molcas":
        mh_en_gr = mh.read_out(in_pos,
                               dyn_bool = True,
                               in_mol = mol_atoms,
                               natoms_flex = natoms_flex, # FJH
                               natoms = natoms,
                               state = state,
                               states = states,
                               mult = mult,
                               singlestate = singlestate,
                               soc_coupling = soc_coupling)
         
#        mh_en_gr = mh.read_out(in_pos)

    # Format energies and gradients for dynamics processing
    """ data format
    mh_en_gr
        0 energy    (nstate, )
        1 grad      (nstate, natoms, 3)
        2 gr_energy float
        3 nac       (nnac, natoms, 3)
        4 soc       (nsoc,)

    ml_en_gr
        0 energy    float
        1 grad      (natoms * 3,)
        2 gr_energy float

    rl_en_gr      
       	0 energy    float
       	1 grad 	    (natoms * 3,)
       	2 gr_energy float
    """
    mh_en, mh_gr, mh_scf, nac, soc = mh_en_gr


    # combine ONIOM result without cap atoms
    # note for future
    # if cap atoms are included, mg_gr and mh_gr need to multiply the jocabian between capped and uncapped model
    en_combo =  mh_en 
    gr_combo =  mh_gr 
    scf_combo = mh_scf
    en_out = en_combo
    gr_out = gr_combo
    nac_out = nac
    soc_out = soc

    evconv = 27.2114 # eV / Hartree

    if stop_traj is not None:
        es_gs_gap = (float(en_combo[state-1]) - float(scf_combo)*evconv
        if es_gs_gap < stop_traj:
           sys.exit("Excited state and ground state energies are too close (DE = {} eV)".format(es_gs_gap))

    # print some updates in the output
    out_file.write("------------------------------\n")
    out_file.write("Iteration: " + str(traj.iter) + "\n")
    out_file.write("Excited state energy: {:>28.8f} eV{:>14.8f} au\n".format(
        float(mh_en[state-1])*evconv, float(mh_en[state-1])))
    out_file.write(
        "Ground state energy: {:>29.8f} eV{:>14.8f} au\n".format(
        float(scf_combo) * evconv, float(scf_combo)))
    out_file.write(" Energy Gap: {:>42.8f} eV{:>14.8f} au\n".format(
            (float(en_combo[state-1]) - float(scf_combo))*evconv,
            float(en_combo[state-1]) - float(scf_combo)))

    out_file.flush()

    # Save all necessary data to restart a trajectory
    if traj.iter % traj.chk_stp == 0:
        mh.save_checkpoint(traj.iter,traj.Eini,traj.E,traj.G,traj.Gp,traj.V,traj.R)
 
    return (en_out, gr_out, nac_out, soc_out)


def atoms_to_fromage(atoms):
    """
    Transform a 3D matrix of atomic Cartesian coordinates into a 1D array
    for use by fromage's sequence() function

    """
    return np.reshape(atoms,(-1,))


###################################################
########### Trajectory Class Definition ###########
###################################################

class Trajectory:
    """
    Defines a Trajectory object which contains as attributes all data necessary for 
        ab initio and nonadiabatic dynamics simulations

     Trajectory class methods allow for the propagation of the trajectory in time 
        according to the velocity-Verlet algorithm, in either the NVE ensemble or 
        NVT ensemble, using a Nose-Hoover thermostat.
    """

####################################################
######## Trajectory Initialization Function ########
####################################################

    def __init__(self,init_dict, mol_atoms, shell_atoms=None):
        """
        Initialize a Trajectory object with settings provided by user and parsed
            by initTrajParams(). Below is a table of possible attributes set by
            the user via fromage.in file. Following that is another table of all 
            other parameters used by Trajectory class methods which are not set by 
            the user, and are initialized internally.

        Parameters
        ----------
        init_dict : Dict
            Dict created by initTrajParams() function containing trajectory settings
            read from user input in fromage.in file

        """
        # Read user input settings from init_dict and reformat if necessary
        self.types = init_dict["types"]
        self.R = init_dict["R"]
        self.M = np.reshape(init_dict["M"],(-1,1)) * 1822.8895
        self.V = init_dict["V"]
        self.state = int(init_dict["State"])
        self.natoms = int(init_dict["natoms"])
        self.mult = list(init_dict["mult"])
        self.states = list(init_dict["states"])
        self.nstates = int(init_dict["nstates"])
        self.stepSize = float(init_dict["stepsize"])
        self.T = float(init_dict["T"])
        self.HopType = init_dict["HopType"].lower()
        self.low = init_dict["low_level"].lower()
        self.high = init_dict["high_level"].lower()
        self.out_file = init_dict["out_file"]
        self.statemult = init_dict["statemult"]
        self.nac_coupling = init_dict["nac_coupling"]
       	self.soc_coupling = init_dict["soc_coupling"]
        self.singlestate = init_dict["singlestate"]
        self.shinfo = 'No surface hopping is performed\n'
        if "natoms_flex" in init_dict.keys():
            self.natoms_flex = int(init_dict["natoms_flex"])
        else:
            self.natoms_flex = 0
        if "temp" in init_dict.keys():
            self.temp = float(init_dict["temp"])
            self.Thermo = True
        else:
            self.Thermo = False
        if "check" in init_dict.keys():
            self.check = init_dict["check"]
        else:
            self.check = "trajectory.chk"
        if "verbose" in init_dict.keys():
            self.verbose = int(init_dict["verbose"])
        else:
            self.verbose = 1
        if "gap" in init_dict.keys():
            self.gap = float(init_dict["gap"])
        else:
            self.gap = 0.5
        if "gapsoc" in init_dict.keys():
            self.gapsoc = float(init_dict["gapsoc"])
        else:
            self.gapsoc = 0.5
        if "nactype" in init_dict.keys():
            self.nactype = str(init_dict["nactype"])
        else:
            self.nactype = 'ktdc'
        if "adjust" in init_dict.keys():
            self.adjust = float(init_dict["adjust"])
        else:
            self.adjust = 1
        if "maxhop" in init_dict.keys():
            self.maxh = float(init_dict["maxhop"])
        else:
            self.maxh = 10
        if "decoherence" in init_dict.keys():
            self.deco = float(init_dict["decoherence"])
        else:
            self.deco = 0.1
        if "substeps" in init_dict.keys():
            self.substep = int(init_dict["substeps"])
        else:
            self.substep = 25
        if "integrate" in init_dict.keys():
            self.integrate = int(init_dict["integrate"])
        else:
            self.integrate = 0
        if "reflect" in init_dict.keys():
            self.reflect = int(init_dict["reflect"])
        else:
            self.reflect = 0
        if "chk_stp" in init_dict.keys():
            self.chk_stp = int(init_dict["chk_stp"])
        else:
            self.chk_stp = 25
        if "alter_orbs" in init_dict.keys():
            self.alter_orbs = str(init_dict["alter_orbs"])
        else:
            self.alter_orbs = 0
        if "e_cons" in init_dict.keys():
            self.e_cons = float(init_dict["e_cons"])
        else:
            self.e_cons = 0.00735 # --> 0.2 eV
        if init_dict["dyn_restart"]:
            self.dyn_restart = init_dict["dyn_restart"]
            self.Gp = init_dict["Gp"]
            self.Gp = np.reshape(self.Gp,(-1,len(self.R),3))
            self.iter = init_dict["curr_step"]
            self.Eini = init_dict["Eini"]
        else:
            self.dyn_restart = 0
            self.Gp = np.zeros((self.nstates,self.natoms,3))
            self.iter = 1
            self.Eini = 0.
        if "nprocs" in init_dict.keys():
            self.nprocs = str(init_dict["nprocs"])
        if "stop_traj" in init_dict.keys():
            self.stop_traj = init_dict["stop_traj"]
#        else:
#            self.nprocs = "1"
 
        # Initialize internal trajectory attributes
        self.Rp = np.zeros_like(self.R)
        self.Rpp = np.zeros_like(self.R)
        self.Vs = np.zeros(4)
        self.E = np.zeros(self.nstates)
        self.Ep = np.zeros_like(self.E)
        self.Epp = np.zeros_like(self.E)
        self.G = np.zeros((self.nstates,self.natoms,3))
        self.Gp = np.zeros_like(self.G)
        self.Gpp = np.zeros_like(self.G)
        self.ncouplings = int(self.nstates*(self.nstates-1)/2)     
        self.N = np.zeros((self.ncouplings,self.natoms,3)) # NAC array
        # current nac is fixe to number of state, the next line can be used when 
        # more flexible defination of nac list is available
        # self.N = np.zeros(0)  # NAC array
        self.S = np.zeros(0)  # SOC array
        self.Np = np.zeros_like(self.N) # previous NAC array
        self.Sp = np.zeros_like(self.S) # previous SOC array
        self.t = 0
        self.Aprev = np.zeros((self.nstates,self.nstates),dtype=complex)
        self.Hprev = np.zeros((self.nstates,self.nstates),dtype=complex)
        self.Dprev = np.zeros((self.nstates,self.nstates),dtype=complex)
        self.Acurr = np.zeros_like((self.Aprev),dtype=complex)
        self.Dcurr = np.zeros_like((self.Dprev),dtype=complex)
        self.Hcurr = np.zeros((self.nstates,self.nstates),dtype=complex)
#       delt is the delta time for the probability in the electronic EOM
#       the 41.34137.. factor is because "delt" is defined in a. u.
        self.delt = self.stepSize * 41.341374575751 / self.substep
        self.Ekin = 0.000
        self.Ekinp = 0.000
        self.Hopped = 0
        self.mol_atoms = mol_atoms
        self.shell_atoms = shell_atoms

###################################################
########## Dynamics Function Definitions ##########
###################################################

    def _thermostat(self):
        """
        Scales the kinetic energy and velocities to maintain constant temperature

        """
        # Check that trajectory is running in NVT ensemble
        if self.Thermo:
            # Calculate scaled velocities and kinetic energy from Nose-Hoover chain
            v,vs,ekin = NoseHoover(self.iter, self.natoms, self.V, self.Ekin, self.Vs, self.temp, self.stepSize)
            # Update trajectory's stored values of velocity and kinetic energy
            self.V    = v
            self.Vs   = vs
            self.Ekin = ekin

    
    def _surfacehop(self):
        """
        Compute surface-hopping probability and change trajectory's state if hopping occurs

        """
        # update previous population, energy matrix, and non-adiabatic coupling matrix
        self.Aprev = np.copy(self.Acurr)
        self.Hprev = np.copy(self.Hcurr)
        self.Dprev = np.copy(self.Dcurr)

        # Determine which surface-hopping algorithm to use based on user input
        if self.HopType == 'fssh':
            hop_method = FSSH
        elif self.HopType == 'gsh':
            hop_method = GSH
        else:
            hop_method = NOSH

        # Compute current population, energy matrix, and non-adiabatic coupling matrix
        at, ht, dt, v, hopped, old_state, state, info = hop_method(self)
    
        # Update trajectory's stored parameters
        self.Acurr = at
        self.Hcurr = ht
        self.Dcurr = dt
        self.V  = v
        self.Hopped = hopped
        self.Old    = old_state
        self.state  = state
        self.shinfo = info
        with open(self.out_file,'a+') as out_file:
            out_file.write(self.shinfo)

    def compute_energy(self):
        """
        Compute the energy and gradients of the current geometry and process
        for trajectory-style formatting
        """

        # Compute energy and gradients using sequence function in fro_run.py
        energy, gradients, NACs, SOCs = dynamics_sequence(self)

        # Convert gradients from Eh/Angstrom to Eh/Bohr
        gradients = gradients * 0.529177


        print("in compute_energy")
        print("gradients")
        print(gradients)
        print("")

        return np.array(energy), gradients, NACs, SOCs


    def _propagate(self):
        """
        Shift computed values from "current" time step to "previous", and
        from "previous" to "previous-previous" to make room for next
        time step. Shift needs to be done before the calculations because
        the calculation will overwrite the "current" values.
        """

        # update previous-preivous and previous coordinates, energies
        #     and gradients
        self.Rpp    = self.Rp.copy()
        self.Rp     = self.R.copy()
        self.Ekinpp = self.Ekinp
        self.Ekinp  = self.Ekin
        self.Epp = self.Ep.copy()
        self.Ep = self.E.copy()
        self.Gpp = self.Gp.copy()
        self.Gp = self.G.copy()
        self.Np = self.N.copy()
        self.Ns = self.S.copy()

        return None


####################################################
############## Main Dynamics Function ##############
####################################################


    def run_dynamics(self):
        """
        """

        if self.dyn_restart:
            self.Eini = self.Eini
            self.V = VerletII(self.iter, self.M, self.G, self.Gp, self.V, self.stepSize, self.state)
            self.Gp = self.G
            self.t = float((self.iter-1) * self.stepSize)

        # Note that when self.iter == 1, VerletI and VerletII will not update geometry and velocity
        # They return the geometry and gradient at t = 0, which is equivalen to do a gradient calculation
        while self.t <= self.T:
            # Record energy, gradient, couplings, geometry, velocity
            self._propagate()

            # Update atomic positions
            self.R = VerletI(self.R, self.V, self.G, self.M, self.stepSize, self.state)

            # Compute the energy and gradients for the new position
            self.E, self.G, NACs, SOCs = self.compute_energy()

            # only check phase when nac is computed
            if self.nactype == "nac":
                self.N = self.check_NACs_phase(NACs)

            # this is needed for surface hopping code, but now leave it as an empty array for future development
            self.S = SOCs

            # Update velocities
            self.V = VerletII(self.iter, self.M, self.G, self.Gp, self.V, self.stepSize, self.state)

            # If using thermostat, scale velocity using Nose-Hoover chain
            self._thermostat()

            # If Trajectory is nonadiabatic, check for surface hopping event and
            #    scale kinetic energy if hop occurs
            self._surfacehop()

            # Calculate kinetic energy
            self.Ekin = np.sum(0.5 * np.rot90(self.M) * np.sum(np.square(self.V),axis=1))
    
            if self.iter == 1:
                self.Eini = float(self.E[self.state - 1]) + float(self.Ekin)
            if self.iter > 1:
                self.check_E_conservation()

            # Print current iteration data to checkpoint file
            self.write_to_chk()
#
#            with open("grad",'a+') as wf:
#                wf.write("Gradient " + str(self.iter) + "\n")
#                for atom in self.G:
#                    np.savetxt(wf, atom)
    
            # Update timestep and iteration number
            self.iter += 1
            self.t += self.stepSize

        return None


    # Check there is energy conservation between current and initial steps according to the selected gap
    # e_cons. If there is no energy conservation, the trajectory is stopped.
    def check_E_conservation(self):
        Etotal_curr = float(self.E[self.state - 1]) + float(self.Ekin)
        if int(self.iter) > 2 and np.abs(Etotal_curr - self.Eini) > self.e_cons:
            print("ENERGY IS NOT CONSERVED IN STEP",self.iter,"DE=",np.abs(Etotal_curr - self.Eini),"a.u")
            print("")
            print("Check the active space for the current step and the dyn_restart file to get the info to restart the trajectory")
            self.write_info_restart_dyn()
            import sys
            sys.exit()

#            else:
#                print("Reordering option activated for this set of orbitals:")
#                print(self.alter_orbs)
#                mh.try_reord_act_spc(self.alter_orbs)
        return
#
    def check_NACs_phase(self,NACs):
        '''
         Compute the phase between two NAC vectors by
         getting the inner product between them. If any
         phase is different, it is corrected to preserve it
        '''

        ncoup = len(NACs)
        nphase = np.zeros(ncoup)
        cossine = np.zeros(ncoup)
        NAC_phase_chk = np.zeros_like(self.N)
        for i in range(ncoup):
            tmp_NAC_prev = self.Np[i,:,:].flatten()
            tmp_NAC_curr = NACs[i,:,:].flatten()
            escalar_prev = np.sqrt(np.dot(tmp_NAC_prev,tmp_NAC_prev))
            escalar_curr = np.sqrt(np.dot(tmp_NAC_curr,tmp_NAC_curr))
            escalar_prev_curr = np.dot(tmp_NAC_prev,tmp_NAC_curr)
#
#           get the phase
            eps=1.E-07
            if (np.abs(escalar_prev*escalar_curr) < eps):
                nphase[i] = 1
                cossine[i] = 2.0 #2 means: cossine not defined
            else:
                cossine[i] = escalar_prev_curr / (escalar_prev*escalar_curr)
                if cossine[i] > 0.0:
                    nphase[i] = 1.
                else:
                    nphase[i] = -1.0
        # phase correction
            NAC_phase_chk[i,:,:] = nphase[i] * NACs[i,:,:]   

        if self.verbose > 1:
            print("NACs phase at step",self.iter,file=open("NACs_phase","a+"))
            print(cossine,file=open("NACs_phase","a+"))
            print(nphase,file=open("NACs_phase","a+"))
                 
        return NAC_phase_chk    
        
            
###################################################
######### Checkpoint File Writing Section #########
###################################################

    def __writeLines(self):
        """
        Writes some lines to the trajectory output file. For styling and readability only.

        """
        # Open the checkpoint file provided by the user
        with open(self.check,"a+") as chk:
            # Write lines for immutable part of the output table
            chk.write("#-------------------------------------------------------------------------")
            # Write lines for each additional column in table
            #    one for each state in the state-averaging
            for state in self.E:
                chk.write("----------------")
            for state in self.E: # Repeated to write populations vs time
                chk.write("----------------")
            chk.write("\n")


    def __writeHeader(self):
        """
        Writes the header line to a trajectory output file. Formats header to include computed
        energies for each electronic state in the state-averaging scheme and the populations of
        each state at every time

        """
        # Open the checkpoint file provided by the user
        with open(self.check,"a+") as chk:
            # Formatting string for immutable part of header line
            header = "{:^6}{:^12}{:^7}{:^16}{:^16}{:^16}"
            chk.write(header.format("# Step","Time","State","E_total","E_kin","E_pot"))
            # Iterate over all states in state-averaging
            for i,state in enumerate(self.E):
                # Print each state's energy to the checkpoint file
                chk.write("{:^16}".format("E_pot_" + str(i)))
            for i,state in enumerate(self.E):
                # Print each state's population to the checkpoint file
                chk.write("{:^16}".format("Pop_st_" + str(i)))
            chk.write("\n")


    def write_to_chk(self):
        """
        Write energy output for each iteration to a file for later processing

        """
        # On the first iteration, overwrite the checkpoint file if it exists
        if self.iter == 1 and os.path.isfile(self.check):
            os.remove(self.check)
        if self.iter == 1 and os.path.exists("velocities_vs_t.chk"):
            os.remove("velocities_vs_t.chk")

        # Open the checkpoint file provided by the user
        with open(self.check,"a+") as chk:
            # On first iteration, write Header to checkpoint file
            if self.iter == 1:
                self.__writeLines()
                self.__writeHeader()
                self.__writeLines()

            # Formatting string for immutable part of current iteration's table entry
            entry = "{step:^6d}{time:^12.4f}{state:^7d}{E_total:^16.6f}{E_kinetic:^16.6f}{E_potential:^16.6f}"
#            entry = "|{step:^6d}|{time:^12.4f}|{state:^7d}|{E_total:^16.6f}|{E_kinetic:^16.6f}|{E_potential:^16.6f}"
            # Format current iteration's data and write to checkpoint file
            chk.write(entry.format(step = int(self.iter), 
                                   time = float(self.t), 
                                   state = int(self.state), 
                                   E_total = (float(self.E[self.state - 1]) + float(self.Ekin)), 
                                   E_kinetic = float(self.Ekin), 
                                   E_potential = float(self.E[self.state - 1]) ))

            # Iterate over all states in state-averaging and write computed potential
            #    energies to checkpoint file
            for Estate in self.E:
                chk.write("{E_pot_state:^16.6f}".format(E_pot_state = float(Estate)))

            Popstates = ([(np.real(x)) for x in np.diag(self.Acurr)])
            for popstate in Popstates:
                chk.write("{Pop_state:^16.6f}".format(Pop_state = float(popstate)))
            chk.write("\n")
         
     
        print("Step=",int(self.iter),file=open("velocities_vs_t.chk","a+"))
        print("time",float(self.t),"fs",file=open("velocities_vs_t.chk","a+"))
        print(self.V,file=open("velocities_vs_t.chk","a+"))

    def write_info_restart_dyn(self):
        """
        """
        if os.path.exists("dyn_restart"):
            os.remove("dyn_restart")
        print("Checkpoint at step",int(self.iter),file=open("dyn_restart","a+"))
        print("----------",file=open("dyn_restart","a+"))
        print("Total Energy at step 1:",self.Eini,file=open("dyn_restart","a+"))
        print("----------",file=open("dyn_restart","a+"))
        print("Molecular Strucure",file=open("dyn_restart","a+"))
        print(self.R,file=open("dyn_restart","a+"))
        print("----------",file=open("dyn_restart","a+"))
        print("Velocities",file=open("dyn_restart","a+"))
        print(self.V,file=open("dyn_restart","a+"))
        print("----------",file=open("dyn_restart","a+"))
        print("Gradients",file=open("dyn_restart","a+"))
        print(self.G,file=open("dyn_restart","a+"))
        print("----------",file=open("dyn_restart","a+"))
        print("Gradients at previous step",file=open("dyn_restart","a+"))
        print(self.Gp,file=open("dyn_restart","a+"))
        print("----------",file=open("dyn_restart","a+"))

        return 
