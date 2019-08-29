#!/usr/bin/env python3
# Adapted with permission from Michael Dommett's code available at
# https://github.com/mdommett/exciton_coupling
# This is an implementation of the diabatisation method by Arago and Troisi in
# Arago, J. & Troisi, A. Dynamics of the Excitonic Coupling in Organic Crystals. Phys. Rev. Lett. 114, 026402 (2015).

import numpy as np
from sys import exit,argv
import argparse

from fromage.utils.exci_coupling import PDA
from fromage.utils.exci_coupling import xyz
from fromage.utils.exci_coupling import read_g09
from fromage.utils.exci_coupling import CATC
from fromage.utils.exci_coupling import diabatize

def main(args):

    au2ev=27.211396132
    oli_states = list(range(1,args.nmerstates+1))


    ############################
    # Point Dipole Approximation
    ############################
    if args.method.upper()=="PDA":
        if len(args.monomerfiles)==2:
            g09_1=args.monomerfiles[0]
            g09_2=args.monomerfiles[1]
            mol_1=read_g09.read_xyz(g09_1)
            mol_2=read_g09.read_xyz(g09_2)

            coords_1=xyz.xyz_to_matrix(mol_1)
            symbols_1=xyz.symbols_from_xyz(mol_1)
            coords_2=xyz.xyz_to_matrix(mol_2)
            symbols_2=xyz.symbols_from_xyz(mol_2)
            COM_1=PDA.centre_of_mass(symbols_1,coords_1)
            COM_2=PDA.centre_of_mass(symbols_2,coords_2)

            TD_1=read_g09.read_TD(g09_1,1)
            TD_2=read_g09.read_TD(g09_2,1)

            PDA_coupling=PDA.PDA_coupling(TD_1,TD_2,COM_1,COM_2)
            if args.units=="au":
                print("PDA coupling: {:.3f} H".format(PDA_coupling))
            else:
                print("PDA coupling: {:.3f} eV".format(PDA_coupling*au2ev))
        else:
            exit("Error! Two monomer files (mf) of G09 output are needed!")
    ############################

    ############################
    # Coulomb ATC method:
    ############################
    elif args.method.upper()=="CATC":
        if len(args.monomerfiles)==2:
            g09_1=args.monomerfiles[0]
            g09_2=args.monomerfiles[1]
            mol_1=read_g09.read_xyz(g09_1)
            mol_2=read_g09.read_xyz(g09_2)
            coords_1=xyz.xyz_to_matrix(mol_1)
            coords_2=xyz.xyz_to_matrix(mol_2)
            NTO_1=read_g09.read_NTO(g09_1,len(mol_1))
            NTO_2=read_g09.read_NTO(g09_2,len(mol_2))

            CATC_coupling=CATC.CATC_coupling(NTO_1,NTO_2,coords_1,coords_2)

            if args.units=="au":
                print("CATC coupling: {:.3f} H".format(CATC_coupling))
            else:
                print("CATC coupling: {:.3f} eV".format(CATC_coupling*au2ev))
        else:
            exit("Error! Two monomer files (mf) of G09 output are needed!")

    ############################
    # dE Method
    ############################
    elif args.method.upper()=="DE":
        if args.nmerfiles is not None:
            if len(args.nmerfiles)==1:
                g09_1=args.nmerfiles[0]
                ES_1=read_g09.read_ES(g09_1,min(oli_states))
                ES_2=read_g09.read_ES(g09_1,max(oli_states))
                dE_coupling=(ES_2-ES_1)/2
                if args.units=="au":
                    print("dE coupling: {:.3f} H".format(dE_coupling))
                else:
                    print("dE coupling: {:.3f} eV".format(dE_coupling*au2ev))
            else:
                exit("Error! For the DE method, one dimer file (-df) should be specified.")
        else:
            exit("Error! For the DE method, a dimer file (-df) must be specified.")
    ############################
    # Diabatization Method
    ############################

    elif args.method.upper()=="DIA":
        Es = [read_g09.read_ES(args.nmerfiles[-1],i) for i in oli_states]

        if args.property.upper()=="TDM":
            TDs_oli = [read_g09.read_TD(args.nmerfiles[0],i) for i in oli_states]
            TDs_oli = np.array(TDs_oli)
            TDs_mono = [read_g09.read_TD(mono,args.monstate) for mono in args.monomerfiles]
            TDs_mono = np.array(TDs_mono)

            oli_props = TDs_oli
            mon_props = TDs_mono

        elif args.property.upper()=="ATC":
            oli_nat = read_g09.read_natoms(args.nmerfiles[0])
            # for each state
            ATC_oli = [read_g09.read_NTO(i,oli_nat) for i in args.nmerfiles]
            ATC_oli = np.array(ATC_oli)
            monomers_natoms = [read_g09.read_natoms(i) for i in args.monomerfiles]
            ATC_mono = [list(read_g09.read_NTO(i,j))*len(args.monomerfiles) for i,j in zip(args.monomerfiles,monomers_natoms)]

            ATC_mono= np.array(ATC_mono)

            oli_props = ATC_oli
            mon_props = ATC_mono

        H=diabatize.diabatize(oli_props,mon_props,Es)

        J=H[0,1]
        if args.units=="au":
            print("Diabatic Hamiltonian H:\n{}\n".format(H))
            print("Diabatic coupling: {:.3f} a.u.".format(J))
        else:
            print("Diabatic Hamiltonian H:\n{}\n".format(H*au2ev))
            print("Diabatic coupling: {:.3f} eV".format(J*au2ev))
    else:
        exit("ERROR: Unrecognised method")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input",help="Input files",type=str,nargs='*')
    parser.add_argument("-m","--method",help="[PDA] Point Dipole Approximation or \
    [CATC] Coulomb Atomic Transition Charges or\
    [dE] Energy Difference or \
    [DIA] Diabatization",required="True")
    parser.add_argument("-p","--property",help="[TDM] Transition Dipole Moments [TDM] \
    or [ATC] Atomic Transition Charges",default="TDM")
    parser.add_argument("-mf","--monomerfiles",help="The log files of the monomer calculations\
    for the diabatization procedure",nargs='*')
    parser.add_argument("-nf","--nmerfiles",help="The log files of the N-mer calculation",nargs='*')
    parser.add_argument("-ms","--monstate",help="Excited state to use for the monomer", default=1,type=int)
    parser.add_argument("-ns","--nmerstates",help="Excited state of N-mer to use",default=2,type=int)
    parser.add_argument("-u","--units",help="Output unit [ev] electronvolts or [au] Hartrees",type=str,default='eV')
    user_input = argv[1:]
    args = parser.parse_args(user_input)

    main(args)
