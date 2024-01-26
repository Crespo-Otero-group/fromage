## velocity verlet for PyRAIMD
## Jingbai Li Feb 13 2020
## Updated by Jordan Cox, May 25, 2021
## Updated by Federico Hernandez Oct 27, 2022

import numpy as np

fs_to_au = 41.341374575751

def NoseHoover(iter: int, natom: int, V, Ekin, Vs, temp, t):
    """
    Calculate velocity scaling factor from Nose-Hoover chain of thermostats (N=2)

    Parameters
    ----------
    iter : int
        Index of the current iteration, beginning from 1
    natom : int
        Number of atoms in high layer
    V : np.array(float)
        XYZ components of velocity for each atom in system as [[Vx1,Vy1,Vz1],[Vx2,Vy2,Vz2]...]
        units of velocity are Bohr/a.u. of time
    Ekin : float
        Kinetic energy of current iteration in hartree (Eh)    
    Vs : List 
        List of parameters used by Nose-Hoover algorithm as [Q1, Q2, V1, V2], initialized by this
        function on the first iteration
    temp : float
        Temperature of the simulation in Kelvin
    t : float
        Step size of the trajectory in femtoseconds

    Return
    ------
    V : np.array(float)
        Updated XYZ components of velocity for each atom, units are Bohr/a.u. of time
    Vs : List
        Updated list of parameters used by Nose-Hoover algorithm as [Q1, Q2, V1, V2]
    Ekin : float
        Updated kinetic energy in hartree (Eh)

    """
    t *= fs_to_au
    kb       = 3.16881*10**-6

    # On first iteration, set up Vs parameters
    if iter == 1:
        freq = 1/(22*fs_to_au) ## 22 fs to au Hz
        Q1 = 3 * natom * temp * kb / freq**2
        Q2 = temp*kb/freq**2
        Vs = [Q1, Q2, 0, 0]
    # On other iterations, compute scale factor and apply to Ekin and V
    else:
        # Read parameters from Vs
        Q1, Q2, V1, V2 = Vs

        # Compute scaling factor 's' using Nose-Hoover algorithm
        G2  = (Q1 * V1**2 - temp * kb) / Q2
        V2 += G2 * t / 4
        V1 *= np.exp(-V2 * t / 8)
        G1  = (2 * Ekin - 3 * natom * temp * kb) / Q1
        V1 += G1 * t / 4
        V1 *= np.exp(-V2 * t / 8)
        s   = np.exp(-V1 * t / 2)

        # Scale kinetic energy by s^2
        Ekin *= s**2
        # Scale velocities by s
        V *= s

        # Update V1 and V2 parameters in Vs using Nose-Hoover algorithm
        V1 *= np.exp(-V2 * t / 8)
        G1  = (2 * Ekin - 3 * natom * temp * kb) / Q1
        V1 += G1 * t / 4
        V1 *= np.exp(-V2 * t / 8)
        G2  = (Q1 * V1**2 - temp * kb) / Q2
        V2 += G2 * t / 4
        Vs[2] = V1
        Vs[3] = V2
        
    return V,Vs,Ekin


def VerletI(R, V, G, M, t, state):
    """
    Update nuclear positions using the velocity-Verlet algorithm

    Parameters
    ----------
    R : np.array(float)
        XYZ coordinates for atoms in high layer, in Angstroms
        coordinates stored as [[X1,Y1,Z1],[X2,Y2,Z2]...]
    V : np.array(float)
        XYZ components of velocity for each atom in R as [[Vx1,Vy1,Vz1],[Vx2,Vy2,Vz2]...]
        units of velocity at Bohr/a.u. of time
    G : np.array(float)
        XYZ components of energy gradient for each atom in R as [[Gx1,Gy1,Gz1],[Gx2,Gy2,Gz2]...]
        units of gradient are Eh/Bohr
    M : np.array(float)
        Mass of each atom in R as 1D array, units of mass are a.u. of mass (NOT amu)
    t : float
        Step size of trajectory in femtoseconds
    state : int
        Root of current electronic state

    Return
    ------
    R : np.array(float)
        Updated XYZ coordinates of atoms in high layer, in Angstroms
        coordinates stored as [[X1,Y1,Z1],[X2,Y2,Z2]...]

    """

    # Skip the first iteration
    if iter == 1:
        return R

    # Convert timestep in femtoseconds to a.u. of time 
    t *= fs_to_au

    # Use gradients for current state only
    G  = G[state-1]

    # Velocity-Verlet equation
#    R += ((V * t) - ((0.5 * G / M) * t**2)) * 0.529177    
    R += (V * t - 0.5 * G / M * t**2) * 0.529177

    return R


def VerletII(iter: int, M, G, G0, V, t, state):
    """
    Update atomic velocities using the velocity-Verlet algorithm

    Paramters
    ---------
    iter : int
        Index of current trajectory iteration, beginning with 1
    M : np.array(float)
        Mass of each atom in R as 1D array, units of mass are a.u. of mass (NOT amu)
    G : np.array(float)
        XYZ components of energy gradient for each atom in R as [[Gx1,Gy1,Gz1],[Gx2,Gy2,Gz2]...]
        units of gradient are Eh/Bohr    
    G0 : np.array(float)
        XYZ components of energy gradient for each atom from previous time step as 
        [[Gx1,Gy1,Gz1],[Gx2,Gy2,Gz2]...], units of gradient are Eh/Bohr
    V : np.array(float)
        XYZ components of velocity for each atom in R as [[Vx1,Vy1,Vz1],[Vx2,Vy2,Vz2]...]
        units of velocity at Bohr/a.u. of time    
    t : float
        Step size of trajectory in femtoseconds
    state : int
        Root of current electronic state    
    
    Return
    ------
    V : np.array(float)
        Updated XYZ components of velocity for each atom in R as [[Vx1,Vy1,Vz1],[Vx2,Vy2,Vz2]...]
        units of velocity at Bohr/a.u. of time
    
    """

    # Skip velocity update for first iteration
    if iter == 1:
        return V

    # Convert timestep from femtoseconds to a.u. of time
    t *= fs_to_au

    # Use gradients for the current state only
    G0 = G0[state-1]
    G = G[state-1]

    # Velocity-Verlet equation, gradient = -1 * force
    V -= 0.5 * (G0 + G) / M * t

    return V

