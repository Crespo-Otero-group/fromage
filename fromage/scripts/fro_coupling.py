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

if __name__=='__main__':
    au2ev=27.211396132
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--method",help="[PDA] Point Dipole Approximation or \
    [CATC] Coulomb Atomic Transition Charges or\
    [dE] Energy Difference or \
    [DIA] Diabatization",required="True")
    parser.add_argument("-p","--property",help="[TDM] Transition Dipole Moments [TDM] \
    or [ATC] Atomic Transition Charges",default="TDM")
    parser.add_argument("-mf","--monomerfiles",help="The log files of the monomer calculations\
    for the diabatization procedure",nargs=2)
    parser.add_argument("-df","--dimerfiles",help="The log files of the dimer calculation",nargs='*')
    parser.add_argument("-ms","--monstate",help="Excited state to use for the monomer", default=1,type=int)
    parser.add_argument("-ds","--dimerstates",help="Excited state of dimer to use",nargs=2,default=[1,2],type=int)
    parser.add_argument("-u","--units",help="Output unit [ev] electronvolts or [au] Hartrees",type=str,required="True")


    parser.add_argument("input",help="Input files",type=str,nargs='*')
    user_input = argv[1:]
    args = parser.parse_args(user_input)

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
        if args.dimerfiles is not None:
            if len(args.dimerfiles)==1:
                g09_1=args.dimerfiles[0]
                ES_1=read_g09.read_ES(g09_1,min(args.dimerstates))
                ES_2=read_g09.read_ES(g09_1,max(args.dimerstates))
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
        if len(args.monomerfiles)==2:
            if args.property.upper()=="TDM":
                dimer=args.dimerfiles[0]
                monomer_A=args.monomerfiles[0]
                monomer_B=args.monomerfiles[1]
                TD_1=read_g09.read_TD(dimer,min(args.dimerstates))
                TD_2=read_g09.read_TD(dimer,max(args.dimerstates))
                E_1=read_g09.read_ES(dimer,min(args.dimerstates))
                E_2=read_g09.read_ES(dimer,max(args.dimerstates))
                TD_A=read_g09.read_TD(monomer_A,args.monstate)
                TD_B=read_g09.read_TD(monomer_B,args.monstate)

                H=diabatize.diabatize(TD_1,TD_2,TD_A,TD_B,E_1,E_2)
                DIA_J=H[0,1]
                if args.units=="au":
                    print("Diabatic coupling: {:.3f} eV".format(DIA_J))
                else:
                    print("Diabatic coupling: {:.3f} eV".format(DIA_J*au2ev))

            elif args.property.upper()=="ATC":
                if len(args.dimerfiles)==3:

                    dimer=args.dimerfiles[0]
                    dimer_state_1=args.dimerfiles[1]
                    dimer_state_2=args.dimerfiles[2]

                    monomer_A=args.monomerfiles[0]
                    monomer_B=args.monomerfiles[1]

                    dimer_natoms=read_g09.read_natoms(dimer)
                    monomer_A_natoms=read_g09.read_natoms(monomer_A)
                    monomer_B_natoms=read_g09.read_natoms(monomer_B)

                    ATC_dimer_state_1=read_g09.read_NTO(dimer_state_1,dimer_natoms)
                    ATC_dimer_state_2=read_g09.read_NTO(dimer_state_2,dimer_natoms)

                    E_1=read_g09.read_ES(dimer,min(args.dimerstates))
                    E_2=read_g09.read_ES(dimer,max(args.dimerstates))

                    ATC_monomer_A=read_g09.read_NTO(monomer_A,monomer_A_natoms)
                    ATC_monomer_B=read_g09.read_NTO(monomer_B,monomer_B_natoms)

                    ATC_monomer_AA=np.concatenate((ATC_monomer_A,ATC_monomer_A),axis=0)
                    ATC_monomer_BB=np.concatenate((ATC_monomer_B,ATC_monomer_B),axis=0)

                    H=diabatize.diabatize(ATC_dimer_state_1,ATC_dimer_state_2,ATC_monomer_AA,ATC_monomer_BB,E_1,E_2)
                    DIA_J=H[0,1]
                    if args.units=="au":
                        print("Diabatic coupling: {:.3f} eV".format(DIA_J))
                    else:
                        print("Diabatic coupling: {:.3f} eV".format(DIA_J*au2ev))

        else:
            exit("Error! Two monomer files (-mf) and one dimer file (-df) must\
            be given!")
