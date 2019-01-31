"""
Tools to carry out the diabatization of the adiabatic Hamiltonian to the diabatic
form, producing the exciton coupling J on the off-diagonal.

The diabatization scheme used herein is proposed by Troisi et. al PRL 114, 026402 (2015)
(https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.026402).
The major detail for the implementation can be found in the supplementary information of
the manuscript.
"""
import numpy as np
def diabatize(dims1,dims2,monA,monB,E1,E2):
    """
    Uses the either the TDMs or ATCs of the s1 and s2 states of the dimer and the s1 state of the two monomers, to
    diabatize the adiabatic Hamiltonian of first two excited states (E1 and E1) to the diabatic
    Hamiltonian, where the off diagonal terms are the couplings J

    Accepts 1xn matrices and state energies as inputs

    Parameters
    ----------
    ATCs1,ATCs2,ATCA,ATCB: 1x3 matrices
    E1,E2: floats of the energy of the s1 and s2 states of the dimer
    Returns
    ----------
    2x2 matrix
    """

    dim=np.concatenate((dims1,dims2)).reshape(2,len(dims1))
    mon=np.concatenate((monA,monB)).reshape(2,len(monA))


    M=np.dot(dim,mon.T)

    U,s,Vt= np.linalg.svd(M)

    C=(np.dot(U,Vt)).transpose()

    E=np.matrix(([E1,0],[0,E2]))

    H=np.dot(np.dot(C,E),C.transpose())

    return H
