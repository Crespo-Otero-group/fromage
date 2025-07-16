#!/bin/python
"""
fro_run-linkatom.py

New fromage implementation for performing ONIOM(QM:QM')-EE calculations where the QM:QM' boundary
cuts through covalent bonds. It is similar to fro_run.py with additional functions for adding link
atoms to the model region, redistributing point charges at the QM:QM' boundary, and a Jacobian for
proper treatment of gradients.

Pre-print: https://doi.org/10.26434/chemrxiv-2025-spfvh

Author: Michael Ingham
"""

import os
import subprocess
from datetime import datetime

from scipy.optimize import minimize
import numpy as np

import matplotlib.pyplot as plt

import fromage.io.read_file as rf
from fromage.utils.atom import Atom

from fromage.io.parse_config_file import bool_cast
from fromage.utils import calc
from fromage.utils import array_operations as ao


from fromage.utils.mol import Mol
import numpy as np
import fromage.io.parse_config_file as pcf


def get_z_index(charge_scheme):
    """get index from z-scheme"""
    import re  # for parsing method type

    # Get Z-scheme
    scheme_index = int(re.search(r"Z(\d+)", charge_scheme, re.IGNORECASE).group(1))
    print("Scheme: Z", scheme_index)

    return scheme_index


def detect_M1_atoms(model, point_charges, z_thresh="2.0"):
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
    Mn_atoms = Mol([])  # point charges to remove

    # get connectivity matrix, and set to real region.
    conn_mat = real_region.get_connectivity_mat()

    # set connectivites for m1 atom
    for m1_atom in m1_atoms:
        # get rid of M1

        if m1_atom not in Mn_atoms:
            Mn_atoms.append(m1_atom)

        # get row of connectivity matrix for this atom
        row_index = real_region.get_index_by_pos(m1_atom)
        m1_conn = conn_mat[row_index]

        # only get atoms from this mol - maybe switch off for covalent systems, might be useful though
        if sys_type == "crystal":
            m1_mol = real_region.select(row_index)
        else:
            m1_mol = point_charges.copy()

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
    charge_corr = tot_mn_char / n_charges

    # Add charge correction
    for pc in out_char:
        if pc not in Mn_charges:
            pc.q += charge_corr
        else:
            pc.q = 0

    return out_char


### Redistributed charge and dipole (RCD) scheme
### Author: Michael Ingham 14/11/2023


def calc_dipole(atom1, atom2):
    """Calculate dipole between two charged atoms"""
    # print("R:", round(atom1.dist(atom2),2) )
    return atom1.dist(atom2) * abs(atom2.q - atom1.q)


def group_M2(real_region, m1_atoms, point_charges, m3=False):
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
            if conn == 1:  # if one bond away
                atom = real_region[ind]
                m2_index = point_charges.get_index_by_pos(atom)
                if m2_index != None:
                    m2_atom = point_charges[m2_index]
                    if (
                        real_region.bonded(m1_atom, m2_atom)
                        and m2_atom in point_charges
                        and m2_atom not in m2_tmp
                        and m2_atom not in m1_atoms
                    ):
                        m2_tmp.append(m2_atom)
            ## adding 24/04 -test, can defo make this better
            elif conn == 2 and m3:
                atom = real_region[ind]
                m3_index = point_charges.get_index_by_pos(atom)

                if m3_index != None:
                    m3_atom = point_charges[m3_index]
                    if (
                        m3_atom in point_charges
                    ):  # and m3_atom not in m2_tmp and m2_atom not in m1_atoms:
                        m2_tmp.append(m3_atom)

        m2_atoms.append(m2_tmp)
    return m2_atoms


def run_new_Z_scheme(
    real_region, point_charges, m1_atoms, z_index, redistribute_charge=True
):
    """Run RCD scheme of Lin and Truhlar"""
    from fromage.utils.atom import Atom

    out_char = point_charges.copy()
    init_char = out_char.get_total_charge()

    print("\n\nInitial total charge: ", init_char)
    print("Number of point charges: ", len(out_char))

    if z_index > 1:
        if z_index == 3:
            m3_bool = True
        else:
            m3_bool = False

        print("m1 atoms:", m1_atoms)

        m2_atoms = group_M2(real_region, m1_atoms, point_charges, m3_bool)

        print("m2 atoms: ", m2_atoms)
        m2_tot = Mol(np.concatenate(m2_atoms))

        m2_tot.write_xyz("m2_vis.xyz")

        for mn_char in m2_tot:
            for charge in out_char:

                if charge.very_close(mn_char, thresh=0.1):
                    out_char.remove(charge)

    for mn_char in m1_atoms:
        for charge in out_char:

            if charge.very_close(mn_char, thresh=0.1):
                out_char.remove(charge)

    if redistribute_charge:
        # apply uniform correction to conserve redistributed charges
        final_char = out_char.get_total_charge()
        charge_corr = (final_char - init_char) / len(out_char)

        for pc in out_char:
            pc.q -= charge_corr

    return out_char


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
        # print("Initial total charge: ", out_char )  #m1_atom.q + sum([atom.q for atom in m2_atom_bonded]))
        n_m2_atoms = len(m2_atom_bonded)
        print("Number of redistribution points: ", n_m2_atoms)
        q0 = m1_atom.q / n_m2_atoms
        print("New value of q0: ", q0)
        if redistribute_dipoles:
            q0 *= 2
            print("Redistributing dipoles... q0_RCD: ", q0)

        for m2 in m2_atom_bonded:
            m2_init = m2.copy()
            out_char.remove(m2)

            q0_char = Atom("H", qIn=q0)
            q0_char.set_pos(0.5 * (m1_atom.get_pos() + m2.get_pos()))
            out_char.append(q0_char)

            if redistribute_dipoles:
                # print("Initial M2 charge: ", m2.q)
                m2.q -= q0 / 2  # print("Redistributing dipoles: m2,k: ", m2.q)

            init_dipole = calc_dipole(m1_atom, m2_init)
            final_dipole = calc_dipole(q0_char, m2)

            # redistribute to q0
            if q0_char not in out_char:
                out_char.append(q0_char)

            # update M2 (if RCD)
            out_char.append(m2)

            print("Initial M1-M2 dipole:", init_dipole)
            print("Final q0-M2 dipole", final_dipole)
            print(
                "change in dipole: ",
                round(100 * (init_dipole - final_dipole) / init_dipole),
                "%",
            )
        print(
            "Final total charge: ",
            out_char.get_total_charge() + sum([atom.q for atom in m2_atom_bonded]),
        )
        # print(out_char)

    fin_char = out_char.get_total_charge()
    print("RC scheme completed")
    print("\n\nFinal total charge: ", fin_char)
    print("Number of point charges: ", len(out_char))

    print("Total change in charge", fin_char - init_char)

    out_char.write_xyz("rc_charge_vis.xyz")
    return out_char


def redistribute_charges(
    region_1,
    in_char,
    real_atoms,
    z_scheme="Z2",
    z_thresh=2.0,
    redistribute_charge=True,
):
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

    # get output file
    here = os.getcwd()
    output_file = open(here + "/prep.out", "a")
    output_file.write(f"\n#### Redistribution scheme: {z_scheme} #####\n")
    output_file.write(
        f"Initial total charge of electronic embedding: {in_char.get_total_charge()}\n"
    )

    # detect problematic charges based of distance
    M1_atoms = detect_M1_atoms(region_1, in_char, z_thresh)
    M1_atoms.write_xyz("M1_atoms.xyz")
    output_file.write(f"Number of M1 atoms detected: {len(M1_atoms)}\n")

    # run charge Zint scheme
    if z_scheme.upper() in ["Z1", "Z2", "Z3"]:
        z_index = int(get_z_index(z_scheme))
        print("M1 atoms: ", M1_atoms)
        out_char = run_new_Z_scheme(
            real_atoms,
            in_char,
            M1_atoms,
            z_index,
            redistribute_charge=redistribute_charge,
        )

    # or RCD
    elif z_scheme.lower() == "rcd":
        out_char = run_RC_scheme(
            real_atoms, in_char, M1_atoms, redistribute_dipoles=True
        )

    elif z_scheme.lower() == "rc":
        out_char = run_RC_scheme(real_atoms, in_char, M1_atoms)

    output_file.write(
        f"Final total charge of electronic embedding: {out_char.get_total_charge()}\n\n"
    )
    output_file.close()

    return out_char


def create_jacobian(n_real_atoms):
    """initialise empty matrix of 3N by 3N"""
    matrix = np.identity(3 * n_real_atoms)
    print("Jacobian shape: ", np.shape(matrix))
    return matrix


def resize_jacobian(matrix_in):
    """increment jacobain; increase each dimension by 3"""
    coord_num = np.shape(matrix_in)[0]
    matrix = np.zeros((coord_num + 3, coord_num + 3))
    matrix[:coord_num, :coord_num] = matrix_in[:coord_num, :coord_num]
    print("Jacobian shape: ", np.shape(matrix))
    return matrix


def derivative(jacobian, con, link, prefactor, xyz1, xyz2, dist2=None, dynamic=False):
    """update jacobian derivatives"""

    # xyz for each con and link
    con3 = con * 3
    link3 = link * 3

    ## counters
    counter1 = np.zeros(3, dtype=int)
    counter2 = np.zeros(3, dtype=int)
    for i in range(3):
        counter1[i] = int(con3 - i + 2)
        counter2[i] = int(link3 - i + 2)
    print("Counter 1: ", counter1)
    print("Counter 2: ", counter2)

    ## nullify counter elemts
    jacobian[counter2[2] :, :] = 0
    jacobian[:, counter2[2] :] = 0

    dyn_sym = 0

    for i in range(3):
        if not dynamic:

            d1 = 1 - prefactor
            d2 = prefactor
            print("Fixed derivative 1: ", d1)
            print("Fixed derivative 2: ", d2)

            jacobian[counter1[i], counter2[i]] = d1
            jacobian[counter2[i], counter2[i]] = d2
            dyn_sym += d2 + d1
        else:

            # calculate derivative with respect to LAC
            d1 = 1 - prefactor * (1 - (xyz2[i] - xyz1[i]) ** 2 / dist2)

            jacobian[counter1[i], counter2[i]] = d1

            if i == 0:
                j1, j2 = 1, 2  # y and z coordinates
            elif i == 1:
                j1, j2 = 0, 2  # x and z coordinates
            else:
                j1, j2 = 0, 1  # x and y coordinates
            print(f"index: {i}, i=j: {j1,j2}")
            print(f"Pre-factor {prefactor}")
            print(f"dist2 {dist2}")
            # print("xyz1", xyz1)
            # print("xyz2", xyz2)
            d2 = (
                prefactor
                * ((xyz2[j1] - xyz1[j1]) ** 2 + (xyz2[j2] - xyz1[j2]) ** 2)
                / dist2
            )

            jacobian[counter2[i], counter2[i]] = d2

            print("Dynamic derivative 1: ", d1)
            print("Dynamic derivative 2: ", d2)
            dyn_sym += d2 + d1

            #
    print("Dynamic sum:", dyn_sym)
    return jacobian


def newcoord(conXyz, hostXyz, jacobian, i, j, real, dynamic=False):
    # C-H
    dist = 1.084

    # link atom generation
    xyz1 = conXyz.get_pos()
    xyz2 = hostXyz.get_pos()

    if conXyz.elem.upper() == "C" and hostXyz.elem.upper() == "C":
        if not dynamic:
            print("Fixed mode active")
            dist2 = 1.528  # AA
            print(f"dist2: {dist2}. Constant value reference: 1.528 A")
        elif dynamic:
            print("Dynamic mode active")

            dist2 = np.linalg.norm(xyz2 - xyz1)

            print("VALUES FOR SCRIPT: ", conXyz, hostXyz)

            print(f"dist2: {dist2}. Constant value reference: 1.528 A")

    ## fixed mode
    prefactor = dist / dist2
    print("Pre-factor: ", prefactor)

    # get link atom
    coord = Atom("H")
    coord.set_pos(xyz1 + (xyz2 - xyz1) * prefactor)

    # LAH-LA difference
    print("LA-LH distance:", np.linalg.norm(coord.get_pos() - xyz1))
    print("LAC-LA distance:", np.linalg.norm(coord.get_pos() - xyz2))

    # update jacobian
    jacobian = derivative(jacobian, i, j, prefactor, xyz1, xyz2, dist2, dynamic=dynamic)

    return coord, jacobian


def jac_transform(jacobian, gradients, vis=False):
    """transform gradients with Jacobian"""
    gradients_init = gradients.reshape(int(len(gradients) / 3), 3)

    print("Gradients initial: ", np.shape(gradients))
    # gradients = gradients.reshape(-1)  # make into 1D array
    gradients = np.matmul(jacobian, gradients)

    gradients_final = gradients.reshape(int(len(gradients) / 3), 3)

    dJac = gradients_final - gradients_init

    print("Change in gradients:\n", dJac)

    print("Final gradients:\n", gradients)

    if vis:
        # Visualization
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))

        cmin = 0
        cmax = max(dJac.max(), gradients.max())

        im1 = axs[0].imshow(
            np.abs(dJac), vmin=cmin, vmax=cmax, aspect="auto", cmap="viridis"
        )
        axs[0].set_title("dJac")
        fig.colorbar(im1, ax=axs[0])

        im2 = axs[1].imshow(
            np.abs(gradients.reshape(int(len(gradients) / 3), 3)),
            vmin=cmin,
            vmax=cmax,
            aspect="auto",
            cmap="viridis",
        )
        axs[1].set_title("Gradients")
        fig.colorbar(im2, ax=axs[1])

        # plt.show()
    return gradients


def read_charge_dat(charge_path, real_mol):

    from fromage.utils.mol import Mol

    with open(charge_path, "r") as f:

        char_mol = Mol([])
        for line in f:
            line = line.split()

            pos = [line[0], line[1], line[2], line[3]]
            map(float, pos)

            new_atom = Atom("H", pos[0], pos[1], pos[2], pos[3])
            char_mol.append(new_atom)

    count = 0

    real_char = real_mol.copy()
    for atom_char in char_mol:
        for atom_real in real_char:
            if atom_char.very_close(atom_real, thresh=0.5):
                count += 1
                # print(count, atom_real, atom_char)
                atom_real.q = atom_char.q

    print(
        "Number of uncharged molecules: ",
        sum([1 for atom in real_char if atom.q == 0.0]),
    )
    return real_char


def get_aug_model(real, model, model_indices):

    aug_model = model.copy()
    jacobian = create_jacobian(len(model))

    nla = 0
    for i in model_indices:
        atom_i = real[i]
        for j in range(len(real)):
            atom_j = real[j]
            if real.bonded(atom_i, atom_j):
                ## add link atom
                if atom_i in model and atom_j not in model:
                    nla += 1
                    print("Adding link atom ", nla, atom_i, atom_j)
                    jacobian = resize_jacobian(jacobian)
                    coord, jacobian = newcoord(atom_i, atom_j, jacobian, i, j, real)

                    print("New link atom: ", coord)
                    aug_model.append(coord)

    print("Number of link atoms: ", nla)
    return aug_model, jacobian


def prep_model(real, model_indices):
    """rearrange model to be at start"""
    from fromage.utils.mol import Mol

    # get model region
    model = Mol([])
    for i in range(max(model_indices) + 1):
        model.append(real[i])
    model.set_bonding(real.bonding, real.thresh)

    # get shell
    shell = real.copy()
    for atom in real:
        for atom_b in model:
            if atom.very_close(atom_b, thresh=0.1):
                shell.remove(atom)

    # make modified real and detect bond cut atoms
    real = model + shell
    lac, lah = shell.detect_bondcuts(model)

    # remove lac in model
    for atom_b in lac:
        for atom in model:
            if atom.very_close(atom_b, thresh=0.5):
                model.remove(atom)

    # move lac to start
    model = Mol([atom for atom in lac]) + model

    # remove lah from shell
    for atom_b in lah:
        for atom in shell:
            if atom.very_close(atom_b, thresh=0.5):
                shell.remove(atom)

    if os.path.exists("flex.xyz"):
        flex = rf.mol_from_file("flex.xyz")
        # remove lah from shell

        for atom in model:
            for atom_b in flex:
                if atom.very_close(atom_b, thresh=0.5):
                    flex.remove(atom)

        for atom in flex:
            for atom_b in shell:
                if atom.very_close(atom_b, thresh=0.5):
                    shell.remove(atom)

        for atom in lah:
            for atom_b in flex:
                if atom.very_close(atom_b, thresh=0.5):
                    flex.remove(atom)

        shell = Mol([atom for atom in lah]) + flex + shell

    else:
        # move lah to start of shell
        shell = Mol([atom for atom in lah]) + shell

    real = model + shell

    model.write_xyz("model_in.xyz")
    real.write_xyz("real_in.xyz")
    shell.write_xyz("shell_in.xyz")

    print("model atoms:", len(model))
    print("real atoms: ", len(real))
    print("shell atoms: ", len(shell))
    print("lac atoms: ", len(lac))
    print("lah atoms:", len(lah))

    num_la = len(lah)

    return model, real, num_la


def singlepoint(atom_array):
    """configure job"""
    global iteration
    iteration += 1

    print("atoms in atom_array: ", len(atom_array) / 3)
    # update atom_array
    atoms_array = np.concatenate([atom_array, end_atoms])
    print("atoms in atoms_array: ", len(atoms_array) / 3)
    # initialise calculation objects
    rl = calc.setup_calc("rl", low_level)
    ml = calc.setup_calc("ml", low_level)
    mh = calc.setup_calc("mh", high_level)

    ## get new real coordinates
    for atom, pos in zip(real, atoms_array.reshape(int(len(atoms_array) / 3), 3)):
        atom.set_pos(pos)

    # update model region
    for atom_a, atom_b in zip(model, real):
        atom_a.set_pos(atom_b.get_pos())

    ## get model and jacobian
    aug_model, jacobian = get_aug_model(real, model, model_indices)  #

    if not os.path.exists("jac.png"):
        plt.imshow(jacobian)

        plt.savefig("jac.png")

    model_array = np.concatenate([atom.get_pos() for atom in aug_model])

    rl_proc = rl.run(ao.array2atom(real, atoms_array), None)
    rl_proc.wait()

    rl_charges = rf.mol_from_file("rl/geom.xyz")

    ## generate point charge embedding between optimisation steps
    recalculate_charge = False
    if not recalculate_charge and iteration == 1:
        with open("fromage.out", "a") as f:
            f.write(
                "No recalculation of PCE. Point charges are fixed to their initial values\n"
            )
        # molden_char = rl.read_charges()
        subprocess.Popen(["cp", f"rl/{charge_keyword}", f"rl/{charge_keyword}_init"])
    elif not recalculate_charge:
        with open("fromage.out", "a") as f:
            f.write("Moving fixed value charges\n")
        subprocess.Popen(["cp", f"rl/{charge_keyword}_init", f"rl/{charge_keyword}"])
        # molden_char = rl.read_charges()
        # print("MOLDEN: ", molden_char)

    subprocess.run(f"head rl/{charge_keyword}", shell=True)
    with open(f"rl/{charge_keyword}", "r") as f:
        molden_char = [float(char) for char in f.readlines()]
        print(molden_char)

    print("charges: ", molden_char)
    with open("molden.char", "a") as f:
        f.write(str(molden_char[:5]) + "\n")

    for atom, char in zip(rl_charges, molden_char):
        # print(atom,char)
        atom.q = float(char)

    model.write_xyz("model_tmp.xyz")

    for atom_b in model:
        for atom_a in rl_charges:
            if atom_b.very_close(atom_a, thresh=0.1):
                rl_charges.remove(atom_a)

    if not scheme == "z0":
        rl_charges = redistribute_charges(
            region_1=model,
            in_char=rl_charges,
            real_atoms=real,
            z_scheme=scheme,
            z_thresh=z_thresh,
            redistribute_charge=True,
        )

    rl_charges.write_xyz("rl_temp.xyz")

    # First, calculate the current charge
    current_charge = sum([atom.q for atom in rl_charges])

    # Write the current charge to the file and read the initial charge
    with open("charge.dat", "a+") as char_file:
        char_file.write(f"{current_charge}\n")
        char_file.seek(
            0
        )  # Move the cursor to the beginning of the file to read the initial charge
        char_init = float(char_file.readline().strip())

    # Constrain total charge to the initial value
    if current_charge != char_init:
        correction = (current_charge - char_init) / len(rl_charges)
        for atom in rl_charges:
            atom.q -= correction

    # Write the corrected charge to a new file
    with open("charge-corrected.dat", "a") as corrected_file:
        corrected_charge = sum([atom.q for atom in rl_charges])
        corrected_file.write(f"{corrected_charge}\n")

    # print("rl_char", rl_char)
    # embedding = ao.array2atom(rl_charges, , [atom.q for atom in rl_charges])
    # [print(atom) for atom in embedding]

    mh_proc = mh.run(ao.array2atom(aug_model, model_array), point_flex=rl_charges)
    mh_proc.wait()

    ml_proc = ml.run(ao.array2atom(aug_model, model_array), point_flex=rl_charges)
    ml_proc.wait()

    ## parse output files
    print("length of model array: ", len(model_array))
    print(model_array)
    ml_en_gr = ml.read_out(model_array)
    mh_en_gr = mh.read_out(model_array)
    rl_en_gr = rl.read_out(atoms_array)
    print("length of grad[1]: ", len(mh_en_gr[1]))
    print(mh_en_gr[1])

    ## ONIOM equation
    en_combo = rl_en_gr[0] - ml_en_gr[0] + mh_en_gr[0]
    scf_combo = rl_en_gr[2] - ml_en_gr[2] + mh_en_gr[2]

    ## ONIOM gradients with jacobian transformation
    ml_en_gr = list(ml_en_gr)
    mh_en_gr = list(mh_en_gr)
    rl_en_gr = list(rl_en_gr)

    rl_en_gr[1] = rl_en_gr[1] * damp_fac
    rl_en_gr = tuple(rl_en_gr)

    if jacobian_transform:

        ml_jac = jac_transform(jacobian, ml_en_gr[1]) * damp_fac
        mh_jac = jac_transform(jacobian, mh_en_gr[1]) * damp_fac

    else:
        ml_jac = ml_en_gr[1]
        mh_jac = mh_en_gr[1]

    ml_en_gr[1] = np.pad(ml_jac, (0, len(rl_en_gr[1]) - len(ml_jac)), "constant")
    mh_en_gr[1] = np.pad(mh_jac, (0, len(rl_en_gr[1]) - len(mh_jac)), "constant")

    ml_en_gr = tuple(ml_en_gr)
    mh_en_gr = tuple(mh_en_gr)

    ## ONIOM gradient equation
    gr_combo = rl_en_gr[1] - ml_en_gr[1] + mh_en_gr[1]

    ## fix some gradients

    # for i, grad in enumerate(gr_combo):
    #     if i > len(aug_model)*3: ## must be 3 for each xyz
    #         print(grad)
    #         gr_combo[i] = 0
    gr_combo = gr_combo[:n_atoms_opt]

    with open(f"{here}/fromage.out", "a") as f:
        f.write("------------------------------\n")
        f.write("Iteration: " + str(iteration) + "\n")
        f.write("Real low energy: {:>30.8f} eV\n".format(rl_en_gr[0] * evconv))
        f.write("Model low energy: {:>29.8f} eV\n".format(ml_en_gr[0] * evconv))
        f.write("Model high energy: {:>28.8f} eV\n".format(mh_en_gr[0] * evconv))
        f.write("ONIOM Total energy: {:>27.8f} eV\n".format(en_combo * evconv))
        f.write("ONIOM SCF energy: {:>29.8f} eV\n".format(scf_combo * evconv))
        f.write(
            "Energy grad. norm: {:>28.8f} eV/A\n".format(
                np.linalg.norm(gr_combo * evconv) / damp_fac
            )
        )
        f.write(
            "Grad RMS: {:>37.8f} eV/A\n".format(
                np.sqrt(np.mean(np.square(gr_combo))) * evconv
            )
        )
        f.write(
            "Max gradient: {:>33.8f} eV/A\n".format(np.max(np.abs(gr_combo)) * evconv)
        )
        f.write("Gap: {:>42.8f} eV\n".format((en_combo - scf_combo) * evconv))
    en_out = en_combo
    gr_out = gr_combo

    with open(f"{here}/geom_out.xyz", "a") as f:
        f.write(f"\n{len(real)}\niteration {iteration}")
        for atom in real:
            f.write(f"\n{atom.elem}   {atom.x}   {atom.y} {atom.z}")

    print(f"{en_out},{gr_out})")
    return (en_out, gr_out)
    # return en_out


if __name__ == "__main__":

    evconv = 27.2114  # Something in Hartree * evconv = Something in eV

    # default settings

    def_inputs = {
        "mol_file": "mol.init.xyz",
        "shell_file": "shell.xyz",
        "out_file": "fromage.out",
        "bool_ci": "0",
        "high_level": "gaussian",
        "jac_off": "0",
        "pyberny": "0",
        "low_level": "gaussian",
        "high_level_mg": None,
        "nprocs": "1",
        "sigma": "3.5",
        "gtol": "1e-5",
        "single_point": "0",
        "dynamics": "0",
        "dyn_restart": "0",
        "relax_qmprime": "0",
        "at_reparam": "0",
        "natoms_flex": "0",
        "bool_la": "0",
        "bonding": "dis",
        "thresh": 1.8,
        "nopt": 0,
        "jac_bool": "1",
        "scheme": "Z3",
        "z_thresh": 1.8,
    }

    ## initialise
    here = os.getcwd()
    inputs = def_inputs.copy()
    if os.path.isfile("fromage.in"):
        new_inputs = rf.read_config("fromage.in")
        inputs.update(new_inputs)

    # useful settings
    bonding = inputs["bonding"]
    thresh = float(inputs["thresh"])
    single_point = bool_cast(inputs["single_point"])
    high_level = inputs["high_level"]
    high_level_mg = inputs["high_level_mg"]
    low_level = inputs["low_level"]
    jac_bool = bool_cast(inputs["jac_bool"])
    z_thresh = float(inputs["z_thresh"])

    # clean up old output
    if os.path.exists("geom_out.xyz"):
        os.remove("geom_out.xyz")
    if os.path.exists("fromage.out"):
        os.remove("fromage.out")

    ## useful input
    jacobian_transform = jac_bool
    nopt = inputs["nopt"]
    scheme = inputs["scheme"]
    if nopt == 0:
        raise ValueError("define number of atoms to optimise")
    else:
        model_indices = [i for i in range(int(nopt))]

    ## write output file
    out_file = "fromage.out"
    with open(out_file, "w") as f:
        start_time = datetime.now()
        f.write("STARTING TIME: " + str(start_time) + "\n")

    ## real file
    real = rf.mol_from_file(inputs["mol_file"])
    real.set_bonding(bonding, thresh)

    # model_indices = [a-1 for a in model_indices]
    model, real, n_la = prep_model(real, model_indices)

    # get coordinates as array
    atoms_array = np.concatenate([atom.get_pos() for atom in real])

    ## number of atoms to be optimised (aug model + lah)

    if os.path.exists("flex.xyz"):
        print("partial relaxation of QM' will be performed")
        flex = rf.mol_from_file("flex.xyz")
        n_atoms_opt = len(flex) * 3

        print("Atoms to optimised: ", n_atoms_opt)

    else:
        n_atoms_opt = (len(model) + n_la) * 3

    print("atoms to be optimised: ", n_atoms_opt / 3)

    atom_array = atoms_array[:n_atoms_opt]
    print("Atoms in atom_array: ", len(atom_array) / 3)
    end_atoms = atoms_array[n_atoms_opt:]
    iteration = 0

    damp_fac = 1  # default to 1

    # char_file = open("molden.char", "w")
    char_file = open("charge.dat", "w")
    # char_file2 = open("charge-corrected.dat", "w")

    if low_level == "xtb_gfnff":
        with open("fromage.out", "a") as f:
            f.write("Running QM/MM with Gfn-FF")
        charge_keyword = "gfnff_charges"
    else:
        charge_keyword = "charges"

    ## ONIOM SCF
    if single_point:
        singlepoint(atom_array)
    else:
        res = minimize(
            singlepoint,
            atom_array,
            jac=True,
            options={"disp": True, "gtol": 1e-5 / damp_fac},
        )  # ,  method="CG")
        # res = minimize(singlepoint, atom_array,options={'disp': True, 'gtol':1e-6}) #, method="Newton-CG")
    with open(out_file, "a") as f:
        f.write("DONE\n")
        end_time = datetime.now()
        f.write("ELAPSED TIME: " + str(end_time - start_time) + "\n")
        f.write("ENDING TIME: " + str(end_time) + "\n")
        f.close()
