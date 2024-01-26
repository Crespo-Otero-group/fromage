"""
This module contains the information from the periodic table
"""

# This list is for reverse searching e.g. having the atomic number but not the
# symbol. The format is:
# Lower case symbol, symbol, atomic number, number of valence electrons,
# covalent radiuys in Angstrom, vdw radius in Angstrom, mass in Da
# Data are taken from the Cambridge Structural Database
periodic_list = [("", "", 0, 0, 0.0, 0.0,0.0),
                 ("point", "point", 0, 0, 0.0, 0.0, 0.0),
                 ("h", "H", 1, 1, 0.23, 1.09, 1.008),
                 ("he", "He", 2, 0, 1.5, 1.4, 4.003),
                 ("c", "C", 6, 4, 0.68, 1.7, 12.011),
                 ("n", "N", 7, 5, 0.68, 1.55, 14.007),
                 ("o", "O", 8, 6, 0.68, 1.52, 15.999),
                 ("f", "F", 9, 7, 0.64, 1.47, 18.998),
                 ("s", "S", 16, 6, 1.02, 1.80, 32.066),
                 ("cu", "Cu", 29, 1, 1.32, 1.4, 63.546),
                 ("zn", "Zn", 30, 2, 1.22, 1.39, 65.39),
                 ("au", "Au", 79, 1, 1.36, 1.66, 196.967)]

periodic = {}

for atom in periodic_list:
    periodic[atom[0]] = {"symbol" : atom[1],
            "at_num" : atom[2],
            "valence_e" : atom[3],
            "cov" : atom[4],
            "vdw" : atom[5],
            "mass" : atom[6]}

bohrconv = 1.88973


def num_to_elem(num):
    """
    Return the element symbol according to input atomic number

    Parameters
    ----------
    num : int
        Atomic number of the element
    """

    for entry in periodic_list:
        if num == entry[2]:
            symbol = periodic[entry[0]]["symbol"]

    return symbol
