"""
This module contains the information from the periodic table
"""

# This list is for reverse searching e.g. having the atomic number but not the
# symbol. The format is:
# Lower case symbol, symbol, atomic number, number of valence electrons,
# vdw radius in Angstrom
periodic_list = [("", "", 0, 0, 0.0),
                 ("point", "point",0, 0, 0),
                 ("h", "H", 1, 1, 1.2),
                 ("c", "C", 6, 4, 1.7),
                 ("n", "N", 7, 5, 1.55),
                 ("o", "O", 8, 6, 1.52),
                 ("f","F",9,7,1.35),
                 ("s","S",16,6,1.80)]

periodic = {}

for atom in periodic_list:
    periodic[atom[0]] = {"symbol": atom[1], "at_num": atom[
        2], "valence_e": atom[3], "vdw": atom[4]}


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
