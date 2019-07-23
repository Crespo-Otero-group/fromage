#!/usr/bin/env python
import numpy as np
import sys
import subprocess
import argparse
from collections import namedtuple as nt

def append_line_to_array(line, array, index):
    """Append line of floats to a numpy array"""
    line_floats = [string_to_float(i) for i in line.split()]
    for num in line_floats:
        array[index] = num
        index += 1
    return index

def string_to_float(in_str):
    """Turn the Fortran format string into a float"""
    if "D" not in in_str:
        in_str = in_str.replace("-","D-")
    out_float = float(in_str.replace("D", "E"))
    return out_float

def read_gauss_rwf(file_name):
    with open(file_name) as file_data:
        reading = False
        ind = 0
        flat_array = None
        for line in file_data:
            if reading:
                ind = append_line_to_array(line, flat_array, ind)
            if "length" in line:
                size = int(line.split()[5])
                reading = True
                flat_array = np.zeros(size)


    return flat_array

def exci_info(file_name):
    """Return information about the excitations"""
    ExciInfo = nt('ExciInfo', 'ncore, nocc_a, nocc_b, nvirt_a, nvirt_b, nex')
    with open(file_name,'r') as in_file:
        for line in in_file:
            if 'NBasis=' in line:
                spl = line.split()
                ncore = int(spl[7])
            if 'NROrb' in line:
                spl = line.split()
                nocc_a = int(spl[3])
                nocc_b = int(spl[5])
                nvirt_a = int(spl[7])
                nvirt_b = int(spl[9])
            if 'roots to seek' in line:
                nex = int(line.split()[6])
                break
    ei = ExciInfo(ncore, nocc_a, nocc_b, nvirt_a, nvirt_b, nex)

    return ei

def trans_coeffs(file_name, exci_info):
    """Return the x and y vectors for alpha and beta

    The order of the coefficients is: 12 lines of garbage, the ket x + ket y
    alpha transitions of the first excitation, then the beta transitions, then
    on to the second excitation and so on. Then the same order (without the 12
    useless lines) for ket x - ket y
    """
    ntrans_a = exci_info.nocc_a * exci_info.nvirt_a
    ntrans_b = exci_info.nocc_b * exci_info.nvirt_b
    ntrans_tot = ntrans_a + ntrans_b

    raw_array = read_gauss_rwf(file_name)
    # list of transition coefficients. Each element is 1 excitation
    exci_coeffs_add = []
    exci_coeffs_sub = []

    for i_ex in range(exci_info.nex):
        a_start = i_ex * ntrans_tot + 12
        a_end = i_ex * ntrans_tot + ntrans_a + 12
        b_start = i_ex * ntrans_tot + ntrans_a + 12
        b_end = i_ex * ntrans_tot + ntrans_tot + 12

        trans_a = raw_array[a_start:a_end]
        trans_b = raw_array[b_start:b_end]

        exci_coeffs_add.append(np.array([trans_a,trans_b]))

    for i_ex in range(exci_info.nex):
        n_skip = int((len(raw_array)-12-exci_info.nex)/2)
        a_start = i_ex * ntrans_tot + 12 + n_skip
        a_end = i_ex * ntrans_tot + ntrans_a + 12 + n_skip
        b_start = i_ex * ntrans_tot + ntrans_a + 12 + n_skip
        b_end = i_ex * ntrans_tot + ntrans_tot + 12 + n_skip

        trans_a = raw_array[a_start:a_end]
        trans_b = raw_array[b_start:b_end]

        exci_coeffs_sub.append(np.array([trans_a,trans_b]))

    ket_x_lis = [(i+j)/2 for i,j in zip(exci_coeffs_add,exci_coeffs_sub)]
    ket_y_lis = [(i-j)/2 for i,j in zip(exci_coeffs_add,exci_coeffs_sub)]


    ket_x_lis = []
    ket_y_lis = []

    for i,j in zip(exci_coeffs_add,exci_coeffs_sub):
        ket_x_lis.append((i[0]+j[0])/2)
        ket_y_lis.append((i[0]-j[0])/2)

    ket_x = np.array(ket_x_lis)
    ket_y = np.array(ket_y_lis)

    return ket_x, ket_y


def main(args):

    # read files
    rwf_file = args.rwf_file
    log_file = args.log_file
    # choose excitation
    excitation_num = args.exci_num - 1

    # gather info from log file
    ei = exci_info(log_file)

    # generate transition coefficients, overlap matrix and mo coefficents files
    subprocess.call("rwfdump "+rwf_file+" trans_coeffs.dat 635R",shell=True)
    subprocess.call("rwfdump "+rwf_file+" overlap.dat 514R",shell=True)
    subprocess.call("rwfdump "+rwf_file+" mo_coeffs.dat 524R",shell=True)

    # get transition coefficients
    x_arr, y_arr = trans_coeffs('trans_coeffs.dat', ei)

    # orbital info
    n_orb_a = ei.nocc_a + ei.nvirt_a
    n_orb_b = ei.nocc_b + ei.nvirt_b
    n_core = ei.ncore
    nbas = n_orb_a + n_core

    # assigning AOs to molecules
    ao_A_labels = list(range(0,int(nbas/2)))
    ao_B_labels = list(range(int(nbas/2),nbas))
    ao_all_labels = list(range(0,nbas))

    laps_array = read_gauss_rwf("overlap.dat")
    # make the overlaps into a lower triangular matrix (with diagonal)
    overlaps = np.zeros((nbas,nbas))
    overlaps[np.tril_indices(nbas)] = laps_array
    # make the matrix square
    overlaps += overlaps.T - np.eye(nbas)

    # get the coefficients for each molecular orbital
    mo_coeffs = read_gauss_rwf("mo_coeffs.dat")
    split_mo_coeffs = np.split(mo_coeffs,n_orb_a+n_core)

    rho_lis = []
    tot_smo = np.zeros((nbas,nbas))
    for mo_coeffs in split_mo_coeffs:
        mo_mat = np.outer(mo_coeffs,mo_coeffs)
        smo_mat = mo_mat * overlaps
        tmp_rho = 0
        for ao_i in ao_A_labels:
            tmp_rho += sum(smo_mat[ao_i])

        rho_lis.append(tmp_rho)
        tot_smo += smo_mat

    # the transitions we care about
    x_trans = x_arr[excitation_num]
    y_trans = y_arr[excitation_num]

    occ_rho = rho_lis[ei.ncore:ei.ncore+ei.nocc_a]
    virt_rho = rho_lis[ei.ncore+ei.nocc_a:]

    sigma_p = 0
    delta_p = 0

    # occupied to virtual transitions
    counter = 0
    for occ_orb in occ_rho:
        for virt_orb in virt_rho:
            tmp_sigma = (virt_orb+occ_orb)*(x_trans[counter]**2)
            tmp_delta = (virt_orb-occ_orb)*(x_trans[counter]**2)
            sigma_p += tmp_sigma
            delta_p += tmp_delta
            counter += 1


    # virtual to occupied transitions
    counter = 0
    for occ_orb in occ_rho:
        for virt_orb in virt_rho:
            tmp_sigma = -(occ_orb+virt_orb)*(y_trans[counter]**2)
            tmp_delta = -(occ_orb-virt_orb)*(y_trans[counter]**2)
            sigma_p += tmp_sigma
            delta_p += tmp_delta
            counter += 1

    sigma_p *= 2
    delta_p *= 2

    classi = ""
    if sigma_p <= 1 - args.thresh:
        classi = "LOC(B)"
    elif 2-sigma_p <= 1 - args.thresh:
        classi = "LOC(A)"
    elif delta_p <= -1 + args.thresh:
        classi = "CT A->B"
    elif 1-delta_p <= args.thresh:
        classi = "CT B->A"
    else:
        classi = "Delocalised"

    print("Excitation: " + str(args.exci_num))
    print("SIGMA A:",sigma_p)
    print("DELTA A:",delta_p)
    print(classi)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("log_file", help="The .log file from the excited state calculation",default="gaussian.rwf")
    parser.add_argument("rwf_file", help="The .rwf file from the excited state calculation",default="gaussian.rwf")
    parser.add_argument("exci_num", help="The excitation to classify",default=1, type=int)
    parser.add_argument("-t", "--thresh", help="Threshold for categorising excitons based on output indices. A higher index means a stricter CT criterion and a lower LOC criterion. Default: 0.5",default=0.5, type=float)
    user_input = sys.argv[1:]
    args = parser.parse_args(user_input)

    main(args)
