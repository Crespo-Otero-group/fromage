"""The class used to manipulate lists of Atoms

This module is split into different files for ease of reading and maintainability.
We define the class in __init__ along with basic functions, which means if the
Mol.__init__() function needs to be invoked by any of the methods, this produces
a circular import. The solution we have chosen is to have in-function imports of
fromage.utils.mol every time this happens.

"""
import numpy as np
from copy import deepcopy

from fromage.utils.atom import Atom
import fromage.io.edit_file as ef


def try_ismol(to_test):
    """ Raise exception if the argument is not a Mol object"""
    if not isinstance(to_test, Mol):
        raise TypeError("Cannot cast " +
                        type(to_test).__name__ + " to Mol object")

def make_mol(atoms):
    """
    Generate a Mol object

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms constituting the Mol
    Returns
    -------
    out_mol : Mol object
        The atoms all in one new Mol object

    """
    out_mol = Mol(atoms)
    return out_mol
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
    vectors : 3 x 3 numpy array
        Lattice vectors of the unit cell
    bonding : string 'dist, 'cov' or 'vdw'
        The method for detecting bonding in this molecule.
        'dis' : distance between atoms < threshold
        'cov' : distance - (cov radius of atom a + of atom b) < threshold
        'vdw' : distance - (vwd radius of atom a + of atom b) < threshold
    thresh : float, optional
        Threshold for the detection. If None, use defaults
    geom_info : GeomInfo object
        Geometry information, including numpy coordinate array, plane coeffs,
        principal and secondary axes.

    """
    from ._listyness import append, extend, insert, remove, index, pop, clear, count, __add__, __len__, __getitem__, __setitem__, __contains__
    from ._bonding import set_bonding, set_bonding_str, bonded, per_bonded
    from ._char import es_pot, change_charges, charges, raw_assign_charges, populate, set_connectivity
    from ._selecting import select, per_select, segregate
    from ._cell_operations import complete_mol, complete_cell, supercell, centered_supercell, trans_from_rad, make_cluster, centered_mols, confined
    from ._geom import GeomInfo, coord_array, calc_coord_array, plane_coeffs, calc_plane_coeffs, axes, calc_axes

    def __init__(self, in_atoms=[], vectors=np.zeros((3, 3)), bonding='dis', thresh=1.8):
        # In case the user feeds a lone atom:
        if isinstance(in_atoms, Atom):
            in_atoms = [in_atoms]
        self.atoms = in_atoms
        self.vectors = vectors
        self.bonding = bonding
        self.thresh = thresh
        self.geom = self.GeomInfo()

    def __repr__(self):
        out_str = ""
        for atom in self.atoms:
            out_str += atom.__str__() + "\n"
        return out_str

    def __str__(self):
        return self.__repr__()

    def copy(self):
        return deepcopy(self)

    def same_atoms_as(self, other_mol):
        """Check if all atoms are the same, even if the oredring differs"""
        result = True
        for atom in self:
            if atom not in other_mol:
                result = False
                break
        return result

    def write_xyz(self, name):
        """Write an xyz file of the Mol"""
        ef.write_xyz(name, self.atoms)

    def empty_mol(self):
        """Return an empty mol with the same properties"""
        new_mol = deepcopy(self)
        new_mol.atoms = []
        return new_mol

    def centroid(self):
        """Return np array of the centroid"""
        N = len(self.atoms)
        centro = np.array([0.0, 0.0, 0.0])
        for atom in self.atoms:
            centro[0] += atom.x
            centro[1] += atom.y
            centro[2] += atom.z
        centro = centro / N
        return centro

    def center_mol(self):
        """Translate molecules to center"""
        cen = self.centroid()
        for atom in self.atoms:
            atom.v_translate(-cen)
        return

    def translate(self, vector):
        """
        Translate Mol by a vector

        Parameters
        ----------
        vector : 3 x 1 numpy array
            Translation vector

        """
        for atom in self.atoms:
            atom.v_translate(vector)
        return

    def translated(self, vector):
        """
        Return translate Mol by a vector

        Parameters
        ----------
        vector : 3 x 1 numpy array
            Translation vector

        Returns
        -------
        out_mol : Mol object
            Translated Mol

        """
        new_mol = self.copy()
        for atom in new_mol.atoms:
            atom.v_translate(vector)
        return new_mol

    def remove_duplicates(self, thresh=0.001):
        """Remove the duplicate atoms"""
        purged_mol = Mol([self.atoms[0]])
        for atom_a in self[1:]:
            unique = True
            for atom_b in purged_mol:
                if atom_a.very_close(atom_b, thresh=thresh):
                    unique = False
                    break
            if unique:
                purged_mol.append(atom_a)
        self.atoms = purged_mol
        return

    def dir_to_frac_pos(self):
        """New mol with atoms in fractional coordinates"""

        out_mol = self.copy()
        # transpose to get the transformation matrix
        M = np.transpose(self.vectors)
        # inverse transformation matrix
        U = np.linalg.inv(M)

        for atom in out_mol:
            # change of basis transformation
            dir_pos = atom.get_pos()
            frac_pos = np.dot(U, dir_pos)
            for i, coord in enumerate(frac_pos):
                # if the coordinate is out of range
                if coord < 0 or coord > 1:
                    # translate it to the range [0,1]
                    frac_pos[i] = coord % 1
            atom.set_pos(frac_pos)
        return out_mol

    def frac_to_dir_pos(self):
        """Move all atoms to direct coordinates"""
        out_mol = self.copy()
        for atom in out_mol:
            new_pos = np.matmul(self.vectors.T, atom.get_pos())
            atom.set_pos(new_pos)

        return out_mol

    def split_in_half(self):
        """
        Split the molecule in half by atom order

        Returns
        -------
        mol_a, mol_b : Mol objects
            Mol objects with the attributes of the original Mol but split in
            half

        """

        if len(self) % 2 != 0:
            raise ValueError(
                "Trying to split a Mol with an odd number of atoms")

        mol_a = self.copy()
        mol_b = self.copy()

        dim_len = len(self)
        mol_a.atoms = self[:int(dim_len / 2)]
        mol_b.atoms = self[-int(dim_len / 2):]

        return mol_a, mol_b
