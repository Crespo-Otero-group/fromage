"""The class used to manipulate lists of Atoms
"""
import numpy as np
# from copy import copy
from copy import deepcopy

from cryspy.utils.atom import Atom
import cryspy.io.edit_file as ef

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
    min_lap : float
        The minimum distance of two overlapping vdw radii for the atoms to be
        considered bonded
    vectors : 3 x 3 numpy array
        Lattice vectors of the unit cell

    """

    def __init__(self, in_atoms, min_lap=0.6, vectors=np.zeros((3, 3))):
        # In case the user feeds a lone atom:
        if isinstance(in_atoms, Atom):
            in_atoms = [in_atoms]
        self.atoms = in_atoms
        self.min_lap = min_lap
        self.vectors = vectors

    def __repr__(self):
        out_str = ""
        for atom in self.atoms:
            out_str += atom.__str__() + "\n"
        return out_str

    def __str__(self):
        return self.__repr__()

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
        return Mol(deepcopy(self).atoms + other_mol.atoms, min_lap = self.min_lap, vectors = self.vectors)

    def __len__(self):
        return len(self.atoms)

    def __eq__(self, other):
        return self.atoms == other.atoms

    def __getitem__(self, index):
        return self.atoms[index]

    def __setitem__(self,index,value):
        self.atoms[index] = value

    def write_xyz(self, name):
        ef.write_xyz(name, self.atoms)

    def select(self, labels):
        """
        Return a molecule out of the current Mol.

        The function returns a new Mol of selected atoms atoms. The selection is
        done by measuring by how much adjacent vdw spheres overlap.

        Parameters
        ----------
        label : int or list of ints
            The number of the atoms from which the molecules are generated.

        Returns
        -------
        selected : Mol object
            The selected molecule
        """

        # Make sure that labels is a list
        if isinstance(labels, int):
            labels = [labels]

        # Check for duplicate labels
        if len(labels) > len(set(labels)):
            raise TypeError("Some labels are repeated")

        selected = Mol(deepcopy([self[i] for i in labels]), min_lap = self.min_lap, vectors = self.vectors)
        remaining = deepcopy(self)
        for atom in selected:
            if atom in remaining:
                remaining.remove(atom)

        old_atoms = deepcopy(selected)
        cont = True
        while cont:
            cont = False
            new_atoms = Mol([])
            for old in old_atoms:
                for rem in remaining:
                    if old.at_lap(rem) >= self.min_lap:
                        new_atoms.append(rem)
                        selected.append(rem)
                        remaining.remove(rem)
                        cont = True
            old_atoms = new_atoms
        return selected
