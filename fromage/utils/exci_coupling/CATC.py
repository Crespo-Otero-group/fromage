"""Contains the tools required to calculate the exciton coupling based on Coulomb
Atomic Transition Charges method
"""
import numpy as np
ang2bohr=1.8897259885789
def CATC_coupling(NTO_1,NTO_2,coordinates_1,coordinates_2):
        """
        Calculates the CATC exciton coupling J based on the Coulomb interaction
        between Atomic Transition Charges in two molecules

        Parameters
        ----------
        NTO_1: List of floats
            List of floats of N-atoms for molecule 1
        NTO_2: List of floats
            List of floats of N-atoms for molecule 2
        coordinates_1: Nx3 array of floats
            Array of x,y,z coordinates for molecule 1
        coordinates_2: Nx3 array of floats
            Array of x,y,z coordinates for molecule 2
        Returns
        ----------
        J: float
            Exciton coupling
        """
        J = 0
        coordinates_1=coordinates_1*ang2bohr
        coordinates_2=coordinates_2*ang2bohr
        for i in range(len(NTO_1)):
            for j in range(len(NTO_2)):
                J+=(NTO_1[i]*NTO_2[j])/np.linalg.norm(coordinates_2[j]-coordinates_1[i])
        return J
