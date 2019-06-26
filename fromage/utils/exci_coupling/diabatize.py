"""
Tools to carry out the diabatization of the adiabatic Hamiltonian to the diabatic
form, producing the exciton coupling J on the off-diagonal.

The diabatization scheme used herein is proposed by Troisi et. al PRL 114, 026402 (2015)
(https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.026402).
The major detail for the implementation can be found in the supplementary information of
the manuscript.
"""
import numpy as np
def diabatize(dimprops,monprops,energies):
    """
    Uses the either the TDMs or ATCs of the s1 and s2 states of the dimer and the s1 state of the two monomers, to
    diabatize the adiabatic Hamiltonian of first two excited states (E1 and E1) to the diabatic
    Hamiltonian, where the off diagonal terms are the couplings J

    Accepts 1xn matrices and state energies as inputs

    Parameters
    ----------
    dimprops : numpy array
        The excited state property for the dimer where dimprops[n] corresponds
        to the nth excited state
    monprops : numpy array
        The excited state property for the monomers where monprops[n] corresponds
        to the nth monomer
    energies : numpy array
        The energies of the dimer in the excited states in increasing order

    Returns
    -------
    H : numpy array
        Diabatic Hamiltonian

    """
    M=np.dot(dimprops,monprops.T)

    U,s,Vt= np.linalg.svd(M)

    C=(np.dot(U,Vt)).transpose()

    E=np.diag(energies)

    H=np.dot(np.dot(C,E),C.transpose())

    return H
