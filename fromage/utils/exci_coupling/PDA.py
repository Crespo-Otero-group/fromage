"""Contains the tools required to calculate the exciton coupling based on Point
Dipole approximation
"""
import numpy as np
from fromage.utils.exci_coupling import elements
def centre_of_mass(symbols,coordinates):
    """
    Calculates the centre of mass (COM) for a set of atomic positions, based on:

    COM_x = m1x1 + m2x2 + ... + mnxn / m1 + m2 + ... + mn
    COM_y = m1y1 + m2y2 + ... + mnyn / m1 + m2 + ... + mn
    COM_z = m1z1 + m2z2 + ... + mnzn / m1 + m2 + ... + mn

    Parameters
    ----------
    symbols: List of strings
        List of elemental symbols
    coordinates: Nx3 array of floats
        Array of x,y,z coordinates
    Returns
    -------
    COM: np.array
        Array of x,y,z component of COM

    """
    if len(symbols)!=len(coordinates):
        exit("Inputs not of the same dimension!")
    masses = np.array([elements.element(i).mass for i in symbols])
    mass_sum = np.sum(masses)

    x_coords=coordinates[:,0]
    y_coords=coordinates[:,1]
    z_coords=coordinates[:,2]
    x_numerator=np.sum(np.dot(masses,x_coords))
    y_numerator=np.sum(np.dot(masses,y_coords))
    z_numerator=np.sum(np.dot(masses,z_coords))

    COM_x = x_numerator/mass_sum
    COM_y = y_numerator/mass_sum
    COM_z = z_numerator/mass_sum
    return np.array([COM_x,COM_y,COM_z])

def PDA_coupling(TD_A,TD_B,COM_A,COM_B):
    """
    Calculates the exciton coupling J based on the following equation:

    J = u1*u2/R12^5 - 3(u1*R12)(R12*u2) / R12^5

    where u1 and u2 are the TD vectors of the two systems and R12 is the vector
    between their centre of masses.

    Parameters
    ----------
    TD_A: np.array
        Array of x,y,z components of transition dipole vector of molecule A
    TD_B: np.array
        Array of x,y,z components of transition dipole vector of molecule B


    coordinates: Nx3 array of floats
        Array of x,y,z coordinates

    Returns
    ----------
    Coupling: Float
        Coupling in atomic units

    """
    ang2bohr=1.8897259885789
    COM_A=COM_A*ang2bohr
    COM_B=COM_B*ang2bohr
    R=np.linalg.norm(COM_B-COM_A)
    return (np.dot(TD_A,TD_B)/R**3) - (3*(np.dot(TD_A,COM_A)*np.dot(TD_B,COM_B))/R**5)
