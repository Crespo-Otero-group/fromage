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
        raise TypeError("Cannot cast " +
                        type(to_test).__name__ + " and Mol object")


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

    def __init__(self, in_atoms, min_lap=0.4, vectors=np.zeros((3, 3))):
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
        return Mol(deepcopy(self).atoms + other_mol.atoms, min_lap=self.min_lap, vectors=self.vectors)

    def __len__(self):
        return len(self.atoms)

    def __eq__(self, other):
        return self.atoms == other.atoms

    def __getitem__(self, index):
        return self.atoms[index]

    def __setitem__(self, index, value):
        self.atoms[index] = value

    def write_xyz(self, name):
        """Write an xyz file of the Mol"""
        ef.write_xyz(name, self.atoms)

    def select(self, labels):
        """
        Return a molecule out of the current Mol.

        The function returns a new Mol of selected atoms atoms. The selection is
        done by measuring by how much adjacent vdw spheres overlap. The returned
        Mol's attributes are new objects obtained via a deep copy.

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

        selected = Mol(deepcopy([self[i] for i in labels]),
                       min_lap=self.min_lap, vectors=self.vectors)
        remaining = deepcopy(self)
        for atom in selected:
            if atom in remaining:
                remaining.remove(atom)

        old_atoms = deepcopy(selected)

        # While there are atoms to add
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
                        cont = True  # An atom was added so continue loop
            old_atoms = new_atoms
        return selected

    def per_select(self, labels, old_pos=False):
        """
        Select a molecule out of a Mol in a periodic system.

        Parameters
        ----------
        labels : int or list of ints
            The number of the atoms from which the molecules are generated
        old_pos : bool
            Option to print the selected molecule at its original coordinates

        Returns
        -------
        selected_img : Mol object
            The atoms belonging to the molecule which is selected with certain
            atoms translated so that the molecule is fully connected without
            periodic boundaries
        selected_old : Mol object (optional)
            The atoms belonging to the molecule which is selected before
            translations

        """

        # Make sure that labels is a list
        if isinstance(labels, int):
            labels = [labels]

        # Check for duplicate labels
        if len(labels) > len(set(labels)):
            raise TypeError("Some labels are repeated")

        # Mol of selected atoms from the unit cell
        selected_old = Mol(deepcopy(
            [self[i] for i in labels]), min_lap=self.min_lap, vectors=self.vectors)
        # Mol of selected atoms where the periodic image
        # atoms are translated back to form a molecule
        selected_img = Mol(deepcopy(
            [self[i] for i in labels]), min_lap=self.min_lap, vectors=self.vectors)

        remaining = deepcopy(self)
        for atom in selected_old:
            if atom in remaining:
                remaining.remove(atom)

        old_atoms = deepcopy(selected_old)

        # While there are atoms to add
        cont = True
        while cont == True:
            cont = False
            new_atoms = Mol([])
            for old in old_atoms:
                for rem in remaining:
                    # contains the distance from the point or image and the
                    # coordinates of the point or image
                    vdw_overlap, per_img = old.per_lap(
                        rem, self.vectors, new_pos=True)
                    # if the atom is close enough to be part of the molecule
                    if vdw_overlap >= self.min_lap:
                        new_atoms.append(per_img)
                        selected_old.append(rem)
                        selected_img.append(per_img)
                        remaining.remove(rem)
                        cont = True  # An atom was added so continue loop
                old_atoms = new_atoms

        if old_pos:
            return selected_img, selected_old
        else:
            return selected_img

    def segregate(self):
        """Separate current Mol in a list of Mols of different molecules"""
        molecules = []  # list of molecules
        remaining = deepcopy(self)

        while len(remaining) > 0:
            molecule = remaining.select(0)
            molecules.append(molecule)
            for atom in molecule:
                remaining.remove(atom)
        return molecules

    def complete_mol(self, labels):
        """
        Take a cell and complete certain molecules

        The objective is to end up with a unit cell where the molecules of interest
        are complete. The rest of the atoms of the cell must remain intact. Note that
        the input atoms are transformed and are the same as are present in the
        output.

        Parameters
        ----------
        labels : int or list of ints
            The number of the atoms from which the molecules are generated
        Returns
        -------
        new_mol : Mol object
            The now complete molecule
        new_cell : Mol object
            The cell with the completed molecule
        """
        new_mol, scattered_mol = self.per_select(labels, old_pos=True)
        new_cell_atoms = deepcopy([a for a in self.atoms if a not in scattered_mol])
        new_cell = deepcopy(self)
        new_cell.atoms = new_cell_atoms

        for atom in new_mol:
            new_cell.append(atom)
        return new_mol, new_cell

    def complete_cell(self):
        """
        Return a cell where atoms have been translated to complete all molecules of
        the cell

        Returns
        -------
        out_cell : Mol object
            The new untruncated cell
        full_mol_l : list of Mol objects
            Each molecule in the untruncated cell

        """
        full_mol_l = []
        remaining = deepcopy(self)

        while len(remaining) != 0:
            full_mol, cell = remaining.complete_mol(0)
            full_mol_l.append(full_mol)
            remaining = cell
            for atom in full_mol:
                if atom in remaining:
                    remaining.remove(atom)

        # Convinently, remaining is now an empty Mol
        out_cell = remaining
        for mol in full_mol_l:
            out_cell.extend(mol)
        return out_cell, full_mol_l
