##### MICHAEL NEW FUNCTIONS FOR Z_SCHEME 13/10/2023
### maybe turn into class

from fromage.utils.mol import Mol
import numpy as np 

def get_z_index(charge_scheme):
    """get index from z-scheme"""
    import re # for parsing method type

    # Get Z-scheme 
    scheme_index = int(re.search(r"Z(\d+)", charge_scheme, re.IGNORECASE).group(1))
    print("Scheme: Z",scheme_index)

    return scheme_index

def detect_M1_atoms(model, point_charges,z_thresh="2.0"):
    """
    Detect potentially problematic point charges within the intermolecular
    bonding threshold, z_thresh
    
    Parameters
    ----------
    model, point_charges : Mol
        model region, and the point charges used to embed it
    z_thresh : float
        intermolecular bonding threshold

    Returns
    -------
    m1_atoms : Mol
        The atoms ass
    
    """
    print(z_thresh)
    m1_atoms = Mol([])
    for atom in model:
        for pc in point_charges:
            distance = atom.dist(pc)
            if distance <= float(z_thresh) and pc not in m1_atoms:
                m1_atoms.append(pc)
    return m1_atoms

def find_Mn_atoms(real_region, point_charges, m1_atoms, z_index, sys_type="crystal"):
    """
    Script which creates Mol object for finding Mn atoms as per the Zn-scheme, 
    where ns the number of bonds to remove. For molecular crystals.

    Paramters
    ---------
    real_region : Mol 
    high_points : Mol
    m1_atoms : Mol
    z_index : float

    Return

    """
    Mn_atoms = Mol([]) # point charges to remove

    #get connectivity matrix, and set to real region.
    conn_mat = real_region.load_connectivity_matrix()

    # set connectivites for m1 atom
    for m1_atom in m1_atoms:
        # get rid of M1 
        
        if m1_atom not in Mn_atoms:
            Mn_atoms.append(m1_atom) 

        # get row of connectivity matrix for this atom
        row_index = real_region.get_index_by_pos(m1_atom)
        m1_conn = conn_mat[row_index]

        #only get atoms from this mol - maybe switch off for covalent systems, might be useful though
        if sys_type=="crystal":
            m1_mol = real_region.select(row_index)
        else:
            m1_mol =point_charges.copy()

        for z_i in range(z_index):
            z_i += 1
            for ind, conn in enumerate(m1_conn): 
                if 0 < conn < z_i:
                    atom = real_region[ind]
                    index = point_charges.get_index_by_pos(atom)
                    if atom not in Mn_atoms and atom in m1_mol:
                        Mn_atoms.append(point_charges[index])

    return Mn_atoms


def run_z_scheme(point_charges, Mn_charges):
    """Redistributes point charges across the cluster according to Zn scheme"""

    # Starting charges
    out_char = point_charges.copy()
    tot_mn_char = Mn_charges.get_total_charge()
    
    # Redistribute Mn charges over rest of cluster
    n_charges = len(out_char) - len(Mn_charges)
    charge_corr = tot_mn_char/n_charges

    #Add charge correction
    for pc in out_char:
        if pc not in Mn_charges:
            pc.q += charge_corr
        else:
            pc.q = 0
    
    return out_char


### Redistributed charge and dipole (RCD) scheme
### Author: Michael Ingham 14/11/2023

def calc_dipole(atom1,atom2):
    """Calculate dipole between two charged atoms"""
    #print("R:", round(atom1.dist(atom2),2) )
    return atom1.dist(atom2)*abs(atom2.q - atom1.q)


def group_M2(real_region, m1_atoms, point_charges):
    """
    Get list of M2 atoms for RC and RCD schemes

    Returns
    -------
    m2_atoms : list of list of Atoms
        M2 atoms grouped in list by which M1 they are bonded to. List in same order 
        as m1_atoms input     
    """
    m2_atoms = []
    conn_mat = real_region.load_connectivity_matrix()
    
    # set connectivites for m1 atom
    for m1_atom in m1_atoms:
        # get row of connectivity matrix for this atom
        row_index = real_region.get_index_by_pos(m1_atom)
        m1_conn = conn_mat[row_index]    
        m2_tmp = []
        for ind, conn in enumerate(m1_conn):    
            if conn == 1: # if one bond away
                atom = real_region[ind]
                m2_index = point_charges.get_index_by_pos(atom)
                if m2_index != None:
                    m2_atom = point_charges[m2_index]
                    if real_region.bonded(m1_atom, m2_atom) and m2_atom in point_charges and m2_atom not in m2_tmp and m2_atom not in m1_atoms:
                        m2_tmp.append(m2_atom)
        m2_atoms.append(m2_tmp)
    return m2_atoms

def run_RC_scheme(real_region, point_charges, m1_atoms, redistribute_dipoles=False):
    """Run RCD scheme of Lin and Truhlar"""
    from fromage.utils.atom import Atom
    out_char = point_charges.copy()
    m2_atoms = group_M2(real_region, m1_atoms, point_charges)

    init_char = out_char.get_total_charge()
    print("\n\nInitial total charge: ", init_char)
    print("Number of point charges: ", len(out_char))

    for m1_atom, m2_atom_bonded in zip(m1_atoms, m2_atoms):
        
        # Elimintate M1 charges
        if m1_atom in out_char:
            out_char.remove(m1_atom)


        print("\n##################")
        print("m1 atom:", m1_atom)
        print("initial m2 atoms:", m2_atom_bonded)
        #print("Initial total charge: ", out_char )  #m1_atom.q + sum([atom.q for atom in m2_atom_bonded]))
        n_m2_atoms = len(m2_atom_bonded)
        print("Number of redistribution points: ", n_m2_atoms)
        q0 = m1_atom.q/n_m2_atoms
        print("New value of q0: ", q0)
        if redistribute_dipoles:
            q0 *= 2
            print("Redistributing dipoles... q0_RCD: ", q0)

        for m2 in m2_atom_bonded:
            m2_init = m2.copy()
            out_char.remove(m2)

            q0_char = Atom("H", qIn=q0)
            q0_char.set_pos(0.5*(m1_atom.get_pos()+m2.get_pos()))
            out_char.append(q0_char)

            if redistribute_dipoles: 
                #print("Initial M2 charge: ", m2.q)
                m2.q -= q0/2                #print("Redistributing dipoles: m2,k: ", m2.q)


            init_dipole = calc_dipole(m1_atom, m2_init)
            final_dipole = calc_dipole(q0_char,m2)
            
            #redistribute to q0
            if q0_char not in out_char:
                out_char.append(q0_char)

            # update M2 (if RCD)
            out_char.append(m2)

            print("Initial M1-M2 dipole:", init_dipole)
            print("Final q0-M2 dipole", final_dipole)
            print("change in dipole: ", round(100*(init_dipole-final_dipole)/init_dipole), "%" )
        print("Final total charge: ",  out_char.get_total_charge() + sum([atom.q for atom in m2_atom_bonded]))
        #print(out_char)
    
    fin_char = out_char.get_total_charge()
    print("RC scheme completed")
    print("\n\nFinal total charge: ", fin_char)
    print("Number of point charges: ", len(out_char))


    print("Total change in charge", fin_char-init_char)

    out_char.write_xyz("rc_charge_vis.xyz")
    return out_char



# def run_RCD_scheme(real_region, point_charges, m1_atoms):
#     """Run RCD scheme of Lin and Truhlar"""
#     from fromage.utils.atom import Atom

#     m2_atoms = []
#     conn_mat = real_region.load_connectivity_matrix()
#     out_char = point_charges.copy()

#     m2_atoms = group_M2(real_region, m1_atoms, point_charges)
#     points_out = Mol([])
    
#     for m1_atom, m2_atom_bonded in zip(m1_atoms, m2_atoms):
#         #remove problematic point charge
#         #out_char.remove(m1_atom)
#         print("\n##################")
#         print("m1 atom:", m1_atom)
#         print("initial m2 atoms:", m2_atom_bonded)
#         print("Initial total charge: ", m1_atom.q + sum([atom.q for atom in m2_atom_bonded]))
#         n_m2_atoms = len(m2_atom_bonded)
#         print("Number of redistribution points: ", n_m2_atoms)
#         q0 = m1_atom.q/n_m2_atoms
#         print("New value of q0: ", q0)
#         q0_rcd = 2*q0
#         print("RCD value of q0: ", q0_rcd)

#         out_char = Mol([])
#         for m2 in m2_atom_bonded:
                
#             q0_char = Atom("H", qIn=q0_rcd)
#             q0_char.set_pos(0.5*(m1_atom.get_pos()+m2.get_pos()))
#             out_char.append(q0_char)

#             m2_rcd = m2.q - q0
#             print("RCD value of M2: ", m2_rcd)
#             print("Change in M2: ", m2.q -m2_rcd)

#             m2_out = m2.copy()
#             m2_out.q = m2_rcd
#             out_char.append(m2_out)

#             init_dipole = calc_dipole(m1_atom, m2)
#             final_dipole = calc_dipole(q0_char,m2_out)
#             print("Initial M1-M2 dipole:", init_dipole)
#             print("Final q0-M2 dipole", final_dipole)
#             print("change in dipole: ", round(100*(init_dipole-final_dipole)/init_dipole), "%" )


#         print("Final total charge: ",  out_char.get_total_charge())
#         print(out_char)
#         for char in out_char:
#             points_out.append(char)

#     points_out.write_xyz("/mnt/c/Users/Michael/Downloads/m1-m2_char.xyz")

#     return point_charges


def redistribute_charges(region_1, in_char, real_atoms, z_scheme="Z2", z_thresh=2.0, sys_type="crystal"):
    """
    Redistribute point charges according to Z-scheme.
    
    Parameters
    ----------
    z_scheme : str
        Number of bonds from Q1 (M1 to M3) up to which point charges should be deleted
        and redistributed. Using more than Z3 is not reccomend
    z_thresh : float
        Intermolecular threshold for identifying M1. Default is 2.0 Angstrom. It should be
        sufficiently large to capture singificant nonbonding interactions which give rise
        to overpolarisation
    in_char : Mol
        point charges for embedding model region
    region_1 : Mol
        model region; used to build real region

    Returns
    -------
    out_char : Mol
        Updated charge objects
    
    """
    import os

    #get output file
    here = os.getcwd()
    output_file = open(here+"/prep.out", "a")
    output_file.write(f"\n#### Redistribution scheme: {z_scheme} #####\n")
    output_file.write(f"Initial total charge of electronic embedding: {in_char.get_total_charge()}\n")    
    
    #detect problematic charges based of distance
    M1_atoms = detect_M1_atoms(region_1,in_char,z_thresh)
    M1_atoms.write_xyz("M1_atoms.xyz")
    output_file.write(f"Number of M1 atoms detected: {len(M1_atoms)}\n")

    # run charge Zint scheme
    if z_scheme in ["z1", "z2", "z3", "Z1", "Z2", "Z3"]:
        z_index = get_z_index(z_scheme)
        Mn_charges = find_Mn_atoms(real_atoms, in_char,M1_atoms,z_index, sys_type=sys_type)
        output_file.write(f"Charges up to M{z_index} to be removed: {len(Mn_charges)}\n")
        #Mn_charges.write_xyz(f"Z{z_index}_charges.xyz")
        #output_file.write(f"Charge to be redistributed: {Mn_atoms.get_total_charge()}\n")
        # delete and redistibuted Mn charges. In principle,  manually custormised point charges could be read in here. 
        out_char = run_z_scheme(in_char, Mn_charges)

    # or RCD 
    elif z_scheme.lower() == "rcd":
        out_char = run_RC_scheme(real_atoms, in_char,M1_atoms, redistribute_dipoles=True)
        
    elif z_scheme.lower() == "rc":
        out_char = run_RC_scheme(real_atoms, in_char,M1_atoms)

    output_file.write(f"Final total charge of electronic embedding: {out_char.get_total_charge()}\n\n")
    output_file.close()

    return out_char

