## Surface hopping for PyRAI2MD
## Python module
## Jingbai Li Feb 22, 2022
## Federico Hernandez Oct, 27, 2022

import sys
import numpy as np

def avoid_singularity(v_i, v_j, i, j):
    ## This fuction avoid singularity of v_i-v_j for i < j 
    ## i < j assumes v_i <= v_j, thus assumes the sign of v_i-v_j is -1

    cutoff = 1e-16

    if i < j:
        sign = -1.0
    else:
        sign = 1.0

    if   v_i == v_j:
        diff = sign*cutoff
    elif v_i != v_j and np.abs(v_i-v_j) < cutoff:
        diff = sign*cutoff
    elif v_i != v_j and np.abs(v_i-v_j) >= cutoff:
        #diff=v_i-v_j
        diff = sign*(v_i-v_j) # force v_i < v_j
    return diff

def Reflect(V, N, reflect):
    ## This function refects velocity when frustrated hopping happens

    if   reflect == 1:
        new_V = -V
    elif reflect == 2:
        new_V = V - 2*np.sum(V*N)/np.sum(N*N)*N

    return new_V

def Adjust(Ea, Eb, V, M, N, adjust, reflect):
    ## This function adjust velocity when surface hopping detected
    ## This function call Reflect if frustrated hopping happens

    Ekin = np.sum(0.5*M*V**2)
    frustrated = 0

    if   adjust == 0:
        dT = Ea-Eb+Ekin
        if dT >= 0:
            f = 1.0
        else:
            new_V = Reflect(V,N,reflect)
            frustrated = 1

    elif adjust == 1:
        dT = Ea-Eb+Ekin
        if dT >= 0:
            f = (dT/Ekin)**0.5
            new_V = f*V
        else:
            new_V = Reflect(V,N,reflect)
            frustrated = 1

    elif adjust == 2:
        a = np.sum(N*N/M)
        b = np.sum(V*N)
        dT = Ea-Eb
        dT = 4*a*dT+b**2
        if dT >= 0:
            if b < 0:
                f = (b+dT**0.5)/(2*a)
            else:
                f = (b-dT**0.5)/(2*a)
            new_V = V - f*N/M
        else:
            new_V = Reflect(V,N,reflect)
            frustrated = 1

    # Added by FJH
    elif adjust == 3:
        a = np.sum(N*N/M) * 0.5
        b = np.sum(V*N)
        dE = Eb-Ea
        Delta = b**2 - 4.*a*dE
        cond1 = np.abs(-b + Delta**0.5)
        cond2 = np.abs(-b - Delta**0.5)
        if Delta >= 0:
            if cond1 < cond2:
                f = (-b + Delta**0.5) / (2.*a)
            else:
                f = (-b - Delta**0.5) / (2.*a)
            V += f*N/M
        else:
            new_V = Reflect(V,N,reflect)
            frustrated = 1
    else:
        new_V = V
    return new_V,frustrated

def dPdt(A, H, D):
    ## This function calculate the gradient of state population
    ## The algorithm is based on Tully's method.John C. Tully, J. Chem. Phys. 93, 1061 (1990)

    ## State density A
    ## Hamiltonian H
    ## Non-adiabatic coupling D, D(i,j) = velo * nac(i,j)
    
    ci = len(A)
    dA = np.zeros((ci, ci),dtype=complex)

    for k in range(ci):
        for j in range(ci):
            for l in range(ci):
                dA[k, j] += A[l, j]*(-1j * H[k, l] - D[k, l]) - A[k, l] * (-1j * H[l, j] - D[l, j])

    return dA

def matB(A, H, D):
    ## This function calculate the B matrix
    ## The algorithm is based on Tully's method.John C. Tully, J. Chem. Phys. 93, 1061 (1990)

    ## State density A
    ## Hamiltonian H
    ## Non-adiabatic coupling D, D(i,j) = velo * nac(i,j)

    ci = len(A)
    b = np.zeros((ci, ci))

    for k in range(ci):
        for j in range(ci):
            b[k, j] = 2 * np.imag(np.conj(A[k, j]) * H[k, j]) - 2*np.real(np.conj(A[k, j]) * D[k, j])

    return b

def kTDC(s1, s2, E, Ep, Epp, dt):
    """ Computing the curvature-driven time-dependent coupling
    The method is based on Truhlar et al J. Chem. Theory Comput. 2022 DOI:10.1021/acs.jctc.1c01080

        Parameters:          Type:
            s1               int        state 1
            s2               int        state 2
            E                ndarray    potential energy in the present step
            Ep               ndarray    potential energy in one step before
            Epp              ndarray    potential energy in two step before
            dt               float      time step

        Return:
            nacme            float      time-dependent nonadiabatic coupling

    """

    dVt = avoid_singularity(E[s1], E[s2], s1, s2)  # s1 < s2, thus dVt <0
    dVt_dt = avoid_singularity(Ep[s1], Ep[s2], s1, s2)  # s1 < s2, thus dVt_dt <0
    dVt_2dt = avoid_singularity(Epp[s1], Epp[s2], s1, s2)  # s1 < s2, thus dVt_2dt <0
    d2Vdt2 = -(dVt - 2 * dVt_dt + dVt_2dt) / dt ** 2  # flip d2Vdt2 to positive
    if d2Vdt2 / dVt > 0:
        nacme = -(d2Vdt2 / dVt) ** 0.5 / 2  # flip nacme to negative due to s1 < s2
    else:
        nacme = 0

    return nacme

def GetNAC(state, new_state, nac_coupling, nac, natom):
    """ Pick up non-adibatic coupling vectors from pre-stored array
        Parameters:          Type:
            state            the current state
            new_state        the new state
            nac_coupling     non-adiabatic coupling pair list
            nac              non-adiabatic coupling array

        Return
            nacv             non-adibatic coupling vectors

    """

    nac_pair = sorted([state - 1, new_state - 1])
    if nac_pair in nac_coupling and len(nac) > 0:
        nac_pos = nac_coupling.index(nac_pair)
        nacv = nac[nac_pos]      # pick up pre-stored non-adiabatic coupling vectors between state and new_state
    else:
        # if the nac vector does not exsit, return an unity matrix
        nacv = np.ones((natom, 3))

    return nacv

def FSSH(traj):
    ## This function integrate the hopping posibility during a time step
    ## This function call dPdt to compute gradient of state population

    A         = traj.Acurr
    H         = traj.Hcurr
    D         = traj.Dcurr
    N         = traj.N
    substep   = traj.substep
    delt      = traj.delt
    iter      = traj.iter
    ci        = traj.nstates
    state     = traj.state
    maxhop    = traj.maxh
    usedeco   = traj.deco
    adjust    = traj.adjust
    reflect   = traj.reflect
    verbose   = traj.verbose
    old_state = traj.state
    new_state = traj.state
    integrate = traj.integrate
    V         = traj.V
    M         = traj.M
    E         = traj.E
    Ekin      = traj.Ekin

    """ Computing the fewest swichest surface hopping
    The algorithm is based on Tully's method.John C. Tully, J. Chem. Phys. 93, 1061 (1990)

        Parameters:          Type:
            traj             class       trajectory class

        Return:              Type:
            At               ndarray     the present state denesity matrix
            Ht               ndarray     the present energy matrix (model Hamiltonian)
            Dt               ndarray     the present nonadiabatic matrix
            Vt               ndarray     the adjusted velocity after surface hopping
            hoped            int         surface hopping decision
            old_state        int         the last state
            state            int         the current(new) state

    """

    A            = traj.Acurr
    H            = traj.Hcurr
    D            = traj.Dcurr
    N            = traj.N
    S            = traj.S
    substep      = traj.substep
    delt         = traj.delt
    iter         = traj.iter
    nstate       = traj.nstates
    state        = traj.state
    maxhop       = traj.maxh
    usedeco      = traj.deco
    adjust       = traj.adjust
    reflect      = traj.reflect
    verbose      = traj.verbose
    old_state    = traj.state
    new_state    = traj.state
    integrate    = traj.integrate
    nactype      = traj.nactype
    V            = traj.V
    M            = traj.M
    E            = traj.E
    Ep           = traj.Ep
    Epp          = traj.Epp
    Ekin         = traj.Ekin
    nac_coupling = traj.nac_coupling
    soc_coupling = traj.soc_coupling
    statemult    = traj.statemult
    At = np.zeros((nstate, nstate), dtype = complex)
    Ht = np.diag(E).astype(complex)
    Dt = np.zeros((nstate, nstate), dtype = complex)
    B = np.zeros((nstate, nstate))
    dB = np.zeros((nstate, nstate))
    dAdt = np.zeros((nstate, nstate), dtype = complex)
    dHdt = np.zeros((nstate, nstate), dtype = complex)
    dDdt = np.zeros((nstate, nstate), dtype = complex)

    hop_g = np.zeros(0)
    hoped = 0
    stop = 0

    ## initialize nac matrix
    if iter > 2:
        for n, pair in enumerate(nac_coupling):
            s1, s2 = pair
            if nactype == 'nac':
                nacme = np.sum(V * N[n]) / avoid_singularity(E[s1], E[s2], s1, s2) 
            elif nactype == 'ktdc':
                nacme = kTDC(s1, s2, E, Ep, Epp, delt * substep)
            Dt[s1, s2] = nacme
            Dt[s2, s1] = -Dt[s1, s2]

    ## initialize soc matrix
    for n, pair in enumerate(soc_coupling):
        s1, s2 = pair
        socme = S[n] / 219474.6  # convert cm-1 to Hartree
        Ht[s1, s2] = socme
        Ht[s2, s1] = socme

    ## initialize state index and order
    stateindex = np.argsort(E)
    stateorder = np.argsort(E).argsort()

    ## start fssh calculation
    if iter < 4:
        At[state - 1, state - 1] = 1
        Vt = V
        info = 'No surface hopping is performed\n'
    else:
        dHdt = (Ht - H) / substep
        dDdt = (Dt - D) / substep
        nhop = 0
        
        if verbose >= 2:
            print('-------------- TEST ----------------')
            print('Iter: %s' % (iter))
            print('Previous Population')
            print(A)
            print('Previous Hamiltonian')
            print(H)
            print('Previous NAC')
            print(D)
            print('Current Hamiltonian')
            print(Ht)
            print('Current NAC')
            print(Dt)
            print('One step population gradient')
            print('dPdt')
            print(dPdt(A,H,D))
            print('matB')
            print(matB(A+dPdt(A,H,D)*delt*substep,H,D)*delt*substep)
            print('Integration start')

        for i in range(substep):
            if integrate == 0:
                B = np.zeros((nstate, nstate))
            g = np.zeros(nstate)
            event = 0
            frustrated=0

            H += dHdt
            D += dDdt

            dAdt = dPdt(A, H, D)
            dAdt *= delt
            A += dAdt
            
            if verbose > 2:
                exceed = []
                for rs in range(nstate):
                    rp = np.diag(np.real(A))[rs]
                    dp = np.abs(np.diag(np.real(dAdt)))[rs]
                    if rp > 1:
                        exceed.append((rp - 1) / (dp + 1e-16))
                    elif rp < 0:
                        exceed.append((0 - rp) / (dp + 1e-16))
                    else:
                        exceed.append(0)

                revert = np.amax(exceed)
                rstate = np.argmax(exceed)

                if revert > 0:
                    print(' numerical instability in state population detected')
                    print(' check A matrix')
                    print(A)
                    print(' check dAdt')
                    print(dAdt)
                    print(' exceed values')
                    print(exceed)
                
            dB = matB(A, H, D)
            B += dB

            for j in range(nstate):
                if j != state - 1:
                    g[j] += np.amax([0, B[j, state - 1] * delt / np.real(A[state - 1, state - 1])])

            z = np.random.uniform(0, 1)

            gsum = 0
            for j in range(nstate):
                gsum += g[stateindex[j]]
                nhop = np.abs(stateindex[j] - state + 1)
                if gsum > z and 0 < nhop <= maxhop:
                    new_state = stateindex[j] + 1
                    event = 1
                    hop_g = np.copy(g)
                    hop_gsum = gsum
                    hop_z = z
                    break

            if verbose > 2:
                print('\nSubIter: %5d' % (i+1))
                print('D nac matrix')
                print(D)
                print('A population matrix')
                print(A)
                print('B transition matrix')
                print(B)
                print('Probabality')
                print(' '.join(['%12.8f' % (x) for x in g]))
                print('Population')
                print(' '.join(['%12.8f' % (np.real(x)) for x in np.diag(A)]))
                print('Random: %s' % (z))
                print('old state/new state: %s / %s' % (state, new_state))

            ## detect frustrated hopping and adjust velocity
            if event == 1:
                NAC = GetNAC(state, new_state, nac_coupling, N, len(V))
                Vt, frustrated = Adjust(E[state - 1], E[new_state - 1], V, M, NAC, adjust, reflect)
                if frustrated == 0:
                    state = new_state

            ## decoherance of the propagation 
            if usedeco != 'OFF':
                deco = float(usedeco)
                tau = np.zeros(nstate)

                ## matrix tau
                for k in range(nstate):
                    if k != state-1:
                        tau[k] = np.abs( 1 / avoid_singularity(
                            np.real(H[state - 1, state - 1]), 
                            np.real(H[k, k]),
                            state - 1,
                            k)) * (1 + deco / Ekin) 

                ## update diagonal of A except for current state
                for k in range(nstate):
                    for j in range(nstate):
                        if k != state - 1 and j != state - 1:
                            A[k, j] *= np.exp(-delt / tau[k]) * np.exp(-delt / tau[j])

                ## update diagonal of A for current state
                Asum = 0.0
                for k in range(nstate):
                    if k != state - 1:
                        Asum += np.real(A[k, k])
                Amm = np.real(A[state - 1, state - 1])
                A[state - 1, state - 1] = 1 - Asum

                ## update off-diagonal of A
                for k in range(nstate):
                    for j in range(nstate):
                        if   k == state - 1 and j != state - 1:
                            A[k, j] *= np.exp(-delt / tau[j]) * (np.real(A[state - 1, state - 1]) / (Amm + 1e-16))**0.5
                        elif k != state - 1 and j == state - 1:
                            A[k, j] *= np.exp(-delt / tau[k]) * (np.real(A[state - 1, state - 1]) / (Amm + 1e-16))**0.5

        ## final decision on velocity
        if state == old_state:   # not hoped
            Vt = V               # revert scaled velocity
            hoped = 0
        else:
            NAC = GetNAC(state, new_state, nac_coupling, N, len(V))
            Vt, frustrated = Adjust(E[old_state - 1], E[state - 1], V, M, NAC, adjust, reflect)

            if frustrated == 0:  # hoped
                hoped = 1
            else:                # frustrated hopping
                hoped = 2
                state = old_state

        At = A

        if len(hop_g) == 0:
            hop_g = g
            hop_z = z
            hop_gsum = gsum

        summary = ''
        for n in range(nstate):
            summary += '    %-5s %-5s %-5s %12.8f\n' % (n + 1, statemult[n], stateorder[n] + 1, hop_g[n])

        info = """Surface hopping information
    Random number:           %12.8f
    Accumulated probability: %12.8f
    state mult  level   probability 
%s
    """ % (hop_z, hop_gsum, summary)

    return At, Ht, Dt, Vt, hoped, old_state, state, info

def GSH(traj):
    """ Computing the fewest swichest surface hopping
        The algorithm is based on Zhu-Nakamura Theory, C. Zhu, Phys. Chem. Chem. Phys., 2014, 16, 25883--25895

        Parameters:          Type:
            traj             class       trajectory class

        Return:              Type:
            At               ndarray     the present state denesity matrix
            Ht               ndarray     the present energy matrix (model Hamiltonian)
            Dt               ndarray     the present nonadiabatic matrix
            Vt               ndarray     the adjusted velocity after surface hopping
            hoped            int         surface hopping decision
            old_state        int         the last state
            state            int         the new state

    """

    iter         = traj.iter
    nstate       = traj.nstates
    state        = traj.state
    verbose      = traj.verbose
    V            = traj.V
    M            = traj.M
    E            = traj.E
    statemult    = traj.statemult
    maxhop       = traj.maxh
    adjust       = traj.adjust
    reflect      = traj.reflect

    # random number
    z = np.random.uniform(0, 1)

    # initialize return values
    old_state = state
    new_state = state
    ic_hop = 0
    is_hop = 0
    hoped = 0
    Vt = V
    hop_type = 'no hopping'

    # initialize state index and order
    stateindex = np.argsort(E)
    stateorder = np.argsort(E).argsort()

    # compute surface hopping probability
    if iter > 2:

        # array of approximate NAC matrix for the same spin multiplicity, unity for different spin
        N = np.ones([nstate, V.shape[0], V.shape[1]])

        # array of hopping probability
        g = np.zeros(nstate)

        # accumulated probability
        gsum = 0

        target_spin = statemult[state - 1]

        for i in range(nstate):

            # skip the present state
            if i == state - 1:
                continue

            state_spin = statemult[i]

            if state_spin == target_spin:
                P, N[i] = InternalConversionProbability(i, traj)
            else:
                P = IntersystemCrossingProbability(i, traj)

            g[i] += P

        event = 0
        for j in range(nstate):
            gsum += g[stateindex[j]]
            nhop = np.abs(stateindex[j] - state + 1)
            if gsum > z and 0 < nhop <= maxhop:
                new_state = stateindex[j] + 1
                event = 1
                break

        # if surface hopping event has occured
        if event == 1:
            # Velocity must be adjusted because hop has occurred
            Vt, frustrated = Adjust(E[old_state - 1], E[new_state - 1], V, M, N[state - 1], adjust = adjust, reflect = reflect)

            # if hop is frustrated, revert the current state to old state
            if frustrated == 1:
                state = old_state
                hoped = 2
            else:
                state = new_state
                hoped = 1

        summary = ''
        for n in range(nstate):
            summary += '    %-5s %-5s %-5s %12.8f\n' % (n + 1, statemult[n], stateorder[n] + 1, g[n])

        info = """Surface hopping information
    Random number:           %12.8f
    Accumulated probability: %12.8f
    state mult  level   probability
%s
    """ % (z, gsum, summary)

    else:
        info = 'No surface hopping is performed\n'

    # allocate zeros vector for population state density
    At = np.zeros([nstate, nstate])

    # assign state density at current state to 1
    At[new_state - 1, new_state - 1] = 1

    # Current energy matrix
    Ht = np.diag(E)

    # Current non-adiabatic matrix
    Dt = np.zeros([nstate, nstate])

    if iter > 2 and verbose >= 2:
        print(info)

    return At, Ht, Dt, Vt, hoped, old_state, new_state, info

def InternalConversionProbability(i, traj):
    """ Computing the probability of internal convertion
        The algorithm is based on Zhu-Nakamura Theory, C. Zhu, Phys. Chem. Chem. Phys., 2014, 16, 25883--25895

        Parameters:          Type:
            i                int         computing state
            traj             class       trajectory class

        Return:              Type:
            P                float       surface hopping probability
            N                ndarray     approximate non-adiabatic coupling vectors

    """

    state        = traj.state
    V            = traj.V
    M            = traj.M
    E            = traj.E
    Ep           = traj.Ep
    Epp          = traj.Epp
    G            = traj.G
    Gp           = traj.Gp
    Gpp          = traj.Gpp
    R            = traj.R
    Rp           = traj.Rp
    Rpp          = traj.Rpp
    Ekinp        = traj.Ekinp
    gap          = traj.gap
    test         = 0

    # determine the energy gap by taking absolute value
    delE = np.abs([E[i] - E[state - 1], Ep[i] - Ep[state - 1], Epp[i] - Epp[state - 1]])

    # total energy in the system at time t2 (t)
    Etotp = Ep[state - 1] + Ekinp

    # average energy in the system over time period
    Ex = (Ep[i] + Ep[state - 1]) / 2

    # early stop if it does not satisfy surface hopping condition
    if np.argmin(delE) != 1 or delE[1] > gap/27.211396132 or Etotp - Ex < 0:
        P = 0
        NAC = np.zeros(V.shape)
        return P, NAC

    dE = delE[1]
    # Implementation of EQ 7
    begin_term = (-1 / (R - Rpp))
    if test == 1: print('IC  EQ 7 R & Rpp: %s %s' % (R, Rpp))
    if test == 1: print('IC  EQ 7 begin term: %s' % (begin_term))
    arg_min = np.argmin([i, state - 1])
    arg_max = np.argmax([i, state - 1])
    if test == 1: print('IC  EQ 7 arg_max/min: %s %s' % (arg_max,arg_min))

    f1_grad_manip_1 = (G[arg_min]) * (Rp - Rpp)
    f1_grad_manip_2 = (Gpp[arg_max]) * (Rp - R)
    if test == 1: print('IC  EQ 7 f1_1/f1_2: %s %s' % (f1_grad_manip_1,f1_grad_manip_1))

    F_ia_1 = begin_term * (f1_grad_manip_1 - f1_grad_manip_2)
    if test == 1: print('IC  EQ 7 done, F_1a_1: %s' % (F_ia_1))

    # Implementation of EQ 8
    f2_grad_manip_1 = (G[arg_max]) * (Rp - Rpp)
    f2_grad_manip_2 = (Gpp[arg_min]) * (Rp - R)
    F_ia_2 = begin_term * (f2_grad_manip_1 - f2_grad_manip_2)
    if test == 1: print('IC  EQ 8 done, F_1a_2: %s' % (F_ia_2))

    # approximate nonadiabatic (vibronic) couplings, which are
    # left out in BO approximation
    NAC = (F_ia_2 - F_ia_1) / (M**0.5)
    NAC = NAC / (np.sum(NAC**2)**0.5)
    if test == 1: print('IC  Approximate NAC done: %s' % (NAC))

    # EQ 4, EQ 5
    # F_A = ((F_ia_2 - F_ia_1) / mu)**0.5
    F_A = np.sum((F_ia_2 - F_ia_1)**2 / M)**0.5
    if test == 1: print('IC  EQ 4 done, F_A: %s' % (F_A))

    # F_B = (abs(F_ia_2 * F_ia_1) / mu**0.5)
    F_B = np.abs(np.sum((F_ia_2 * F_ia_1) / M))**0.5
    if test == 1: print('IC  EQ 5 done, F_B: %s' % (F_B))

    # compute a**2 and b**2 from EQ 1 and EQ 2
    # ---- note: dE = 2Vx AND h_bar**2 = 1 in Hartree atomic unit
    a_squared = (F_A * F_B) / (2 * dE**3)
    b_squared = (Etotp - Ex) * (F_A / (F_B * dE))
    if test == 1: print('IC  EQ 1 & 2 done, a^2, b^2: %s %s' % (a_squared,b_squared))

    # GOAL: determine sign in denominator of improved Landau Zener formula for switching
    # probability valid up to the nonadiabtic transition region
    #F_1 = E[i] - Epp[state - 1] # approximate slopes
    #F_2 = E[state - 1] - Epp[i] # here

    #if (F_1 == F_2):
    #    sign = 1
    #else:
    #    # we know the sign of the slope will be negative if either F_1 or
    #    # F_2 is negative but not the other positive if both positive or both negative
    #    sign = np.sign(F_1 * F_2)
    sign = np.sign(np.sum(F_ia_1 * F_ia_2))
    if test == 1: print('IC  Compute F sign done: %s' % (sign))

    # sign of slope determines computation of surface
    # hopping probability P (eq 3)
    pi_over_four_term = -(np.pi/ (4 *(a_squared)**0.5))
    if test == 1: print('IC  P numerator done: %s' % (pi_over_four_term))
    b_in_denom_term = (2 / (b_squared + (np.abs(b_squared**2 + sign))**0.5))
    if test == 1: print('IC  P denomerator done: %s' % (b_in_denom_term))
    P = np.exp(pi_over_four_term * b_in_denom_term**0.5)
    if test == 1: print('IC  P done: %s' % (P))

    return P, NAC

def IntersystemCrossingProbability(i, traj):
    """ Computing the probability of intersystem crossing
        The algorithm is based on Zhu-Nakamura Theory, C. Zhu, Phys. Chem. Chem. Phys., 2020,22, 11440-11451
        The equations are adapted from C. Zhu, Phys. Chem. Chem. Phys., 2014, 16, 25883--25895

        Parameters:          Type:
            i                int         computing state
            traj             class       trajectory class

        Return:              Type:
            P                float       surface hopping probability

    """

    state        = traj.state
    soc_coupling = traj.soc_coupling
    soc          = traj.Sp
    M            = traj.M
    E            = traj.E
    Ep           = traj.Ep
    Epp          = traj.Epp
    Gp           = traj.Gp
    Ekinp        = traj.Ekinp
    gap          = traj.gapsoc
    test         = 0

    # determine the energy gap and type of crossing
    delE = [E[i] - E[state - 1], Ep[i] - Ep[state - 1], Epp[i] - Epp[state - 1]]
    #parallel = np.sign(delE[0]* delE[2])
    parallel = -1 # assume non-parallel PESs

    # total energy in the system at time t2 (t)
    Etotp = Ep[state - 1] + Ekinp

    # set hopping point energy to target state
    Ex = Ep[i]

    # early stop if it does not satisfy surface hopping condition
    if np.argmin(np.abs(delE)) != 1 or np.abs(delE[1]) > gap/27.211396132 or Etotp - Ex < 0:
        P = 0
        return P

    # early stop if it soc was not computed (ignored)
    soc_pair = sorted([state - 1, i])
    if soc_pair not in soc_coupling:
        P = 0
        return P

    # get soc coupling
    soc_pos = soc_coupling.index(soc_pair)
    if len(soc) >= soc_pos + 1:
        soclength = soc[soc_pos]
    else:
        sys.exit('\n  DataNotFoundError\n  PyRAI2MD: looking for spin-orbit coupling between %s and %s' % (state, i + 1))

    V12x2 = 2 * soclength / 219474.6  # convert cm-1 to hartree

    # Implementation of EQ 7
    F_ia_1 = Gp[state - 1]
    F_ia_2 = Gp[i]
    if test == 1: print('ISC EQ 7 done: %s' % (F_ia_1))
    if test == 1: print('ISC EQ 8 done: %s' % (F_ia_2))

    # EQ 4, EQ 5
    F_A = np.sum((F_ia_2 - F_ia_1)**2 / M)**0.5
    if test == 1: print('ISC EQ 4 done, F_A: %s' % (F_A))

    F_B = np.abs(np.sum((F_ia_2 * F_ia_1) / M))**0.5
    if test == 1: print('ISC EQ 5 done, F_B: %s' % (F_B))

    # compute a**2 and b**2 from EQ 1 and EQ 2
    # ---- note: V12x2 = 2 * SOC AND h_bar**2 = 1 in Hartree atomic unit
    a_squared = (F_A * F_B) / (2 * V12x2**3)
    b_squared = (Etotp - Ex) * (F_A / (F_B * V12x2))
    if test == 1: print('ISC EQ 1 & 2 done: %s %s' % (a_squared,b_squared))

    # GOAL: determine sign in denominator of improved Landau Zener formula for switching
    # probability at corssing region
    sign = np.sign(np.sum(F_ia_1 * F_ia_2))
    if test == 1: print('ISC Compute F sign done: %s' % (sign))

    # hopping probability P (eq 3)
    pi_over_four_term = -(np.pi/ (4 *(a_squared)**0.5))
    if test == 1: print('LZ-P numerator done: %s' % (pi_over_four_term))
    b_in_denom_term = (2 / (b_squared + (np.abs(b_squared**2 + sign))**0.5))
    if test == 1: print('LZ-P denomerator done: %s' % (b_in_denom_term))
    P = np.exp(pi_over_four_term * b_in_denom_term**0.5)
    if test == 1: print('LZ-P done: %s' % (P))
    if test == 1: print("""parallel crossing: %s
 1 - P / (P + 1) = %s
 1 - P           = %s
""" % (parallel, 1 - P / (P + 1), 1 - P))

    if parallel == 1:
        P = 1 - P / (P + 1)
    else:
        P = 1 - P

    return P

def NOSH(traj):
    """
    Fake surface hopping method to do single state molecular dynamics

    """
    ci        = traj.nstates
    state     = traj.state
    old_state = traj.state
    hoped     = 0
    V         = traj.V
    E         = traj.E

    # allocate zeros vector for population state density
    At = np.zeros([ci,ci])

    # assign state density at current state to 1
    At[state - 1, state - 1] = 1

    # Current energy matrix
    Ht = np.diag(E)

    # Current non-adiabatic matrix
    Dt = np.zeros([ci, ci])

    # Return the same velocity
    Vt = V

    return At, Ht, Dt, Vt, hoped, old_state, state

