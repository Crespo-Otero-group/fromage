"""The class used to manipulate lists of Atoms
"""
import numpy as np
from copy import copy

from cryspy.utils.atom import Atom

def try_ismol(to_test):
    """ Raise exception if the argument is not a Mol object"""
    if not isinstance(to_test, Mol):
        raise TypeError("Cannot cast " + type(to_test).__name__ + " and Mol object")

class Mol(object):
    """
    Object representing a list of atoms.

    This class can be used to reflect any number of molecules, point charges or
    unit cells. Although Mol shares many methods with list, it deliberately does
    not inherit it in order to avoid nonsensical operations such as Mol1 > Mol2

    Attributes
    ----------
    atoms : list of Atom objects
        Member atoms of Mol
    max_bl : float
        The maximum distance which is considered a bond for Hydrogen
    vectors : 3 x 3 numpy array
        Lattice vectors of the unit cell

    """

    def __init__(self, in_atoms, max_bl=1.7, vectors=np.zeros((3, 3))):
        # In case the user feeds a lone atom:
        if isinstance(in_atoms, Atom):
            in_atoms = [in_atoms]
        self.atoms = in_atoms
        self.max_bl = max_bl
        self.vectors = vectors

        # For iterating
        self.count = 0
        self.limit = len(self.atoms) - 1

    def __iter__(self):
        return self

    def __next__(self):
        self.count += 1
        if self.count > self.limit:
            raise StopIteration
        else:
            return self.atoms[self.count]

    def __repr__(self):
        out_str = ""
        for atom in self.atoms:
            out_str += atom.__str__() + "\n"
        return out_str

    # list-y behaviour
    def append(self, element):
        self.atoms.append(element)

    def extend(self, other_mol):
        self.atoms.extend(other_mol.atoms)

    def insert(self, i, element):
        self.atoms.insert(i, element)

    def remove(self, element):
        self.atoms.remove(element)

    def pop(self, i=-1):
        return self.atoms.pop(i)

    def clear(self):
        self.atoms.clear()

    def count(self, element):
        return self.atoms.count()

    def __add__(self, other_mol):
        try_ismol(other_mol)
        return (Mol(copy(self).atoms + other_mol.atoms))

    def __len__(self):
        return len(self.atoms)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return self.atoms == other.atoms

    # def select(self, label):
    #     """
    #     Return a molecule out of the current Mol.
    #
    #     The function puts a labelled atom in a list of selected atoms. It then
    #     checks for nearby atoms to the selected atoms. If any are found and are not
    #     already selected, they are added to the list. The loop ends when it is
    #     executed and not atoms are added.
    #
    #     Parameters
    #     ----------
    #     label : int
    #         The number of the atom from which the molecule is generated. Note that
    #         it is the Python label, so starting from 0, unlike the label in most
    #         chemistry software which counts from 1
    #
    #     Returns
    #     -------
    #     selected : list of Atom objects
    #         The atoms belonging to the molecule which is selected
    #
    #     """
    #
    #     selected = [copy(self.atoms[label])]
    #     n = True
    #     while n == True:
    #         n = False
    #         old_atoms = new_atoms
    #         new_atoms = []
    #         for i in old_atoms:
    #             for j in self.atoms:
    #                 # find atoms bonded to the newest atoms in the list
    #                 if i.dist(j.x, j.y, j.z) <= max_r and j not in selected:
    #                     new_atoms.append(j)
    #                     selected.append(j)
    #                     # if no more atoms are found stop the loop
    #                     n = True
    #
    #     return selected

    next = __next__  # For Python 2.X
