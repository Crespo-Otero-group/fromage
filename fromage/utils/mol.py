"""The class used to manipulate lists of Atoms
"""
import numpy as np
# from copy import copy
from copy import deepcopy

from fromage.utils.atom import Atom
import fromage.io.edit_file as ef

def try_ismol(to_test):
    """ Raise exception if the argument is not a Mol object"""
    if not isinstance(to_test, Mol):
        raise TypeError("Cannot cast " +
                        type(to_test).__name__ + " to Mol object")

default_thresh = {'dis' : 1.8,
                'cov' : 0.2,
                'vdw' : -0.3}

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

    """



    def __init__(self, in_atoms=[], vectors=np.zeros((3, 3)), bonding = 'dis', thresh = 1.8):
        # In case the user feeds a lone atom:
        if isinstance(in_atoms, Atom):
            in_atoms = [in_atoms]
        self.atoms = in_atoms
        self.vectors = vectors
        self.bonding = bonding
        self.thresh = thresh
    def __repr__(self):
        out_str = ""
        for atom in self.atoms:
            out_str += atom.__str__() + "\n"
        return out_str

    def __str__(self):
        return self.__repr__()

    def set_bonding(self, bonding='dis', thresh=None):
        """
        Set the type of bonding detection used in this Mol

        Parameters
        ----------
        bonding : string 'dis', 'cov' or 'vdw'
            The method for detecting bonding in this molecule.
            'dis' : distance between atoms < threshold
            'cov' : distance - (cov radius of atom a + of atom b) < threshold
            'vdw' : distance - (vwd radius of atom a + of atom b) < threshold
        thresh : float, optional
            Threshold for the detection. If None, use defaults:
            'dis' -> 1.8
            'cov' -> 0.2
            'vdw' -> -0.3

        """
        if bonding not in default_thresh:
            raise TypeError("Unrecognised bonding type: "+ bonding)
        self.bonding = bonding
        if thresh:
            self.thresh = thresh
        else:
            self.thresh = default_thresh[bonding]
        return

    def set_bonding_str(self, in_str):
        """
        Set the type of bonding and threshold with one string

        The string is of the type "cov0.2" or "dis1.7" etc. But giving just the
        threshold or just the bonding gives the default for the ommitted part.
        The order of the bonding and threshold does not matter, so "vdw2.2" is
        the same as "2.2vdw"

        Parameters
        ----------
        in_str : str
            The string which determines the bonding where the threshold and the
            distance are set to default if none are supplied

        """
        bondings = default_thresh.keys()
        bonding = ''
        thresh = None
        # check if bonding has been specified
        for i_bonding in bondings:
            if i_bonding in in_str:
                bonding = i_bonding
        # if there is bonding, try to find threshold
        if bonding:
            stripped = in_str.replace(bonding,'')
            # if there is still a thresh in this string
            if stripped:
                thresh = float(stripped)
        # if only the thresh is specified
        else:
            thresh = float(in_str)
        # if both present
        if bonding and thresh:
            self.set_bonding(bonding=bonding,thresh=thresh)
        # if only bonding
        if bonding and not thresh:
            self.set_bonding(bonding=bonding)
        # if only thresh
        if thresh and not bonding:
            self.set_bonding(thresh=thresh)
        return

    def bonded(self, atom_a, atom_b):
        """
        Check if atom_a is bonded to atom_b given the bonding settings

        Parameters
        ----------
        atom_a, atom_b : Atom objects
            The atoms to be compared
        Returns
        -------
        bonded_bool : bool
            True if the atoms are bonded and False if not
        """
        bonded_bool = atom_a.dist(atom_b, ref=self.bonding) <= self.thresh
        return bonded_bool

    def per_bonded(self, atom_a, atom_b):
        """
        Check if atom_a is bonded to atom_b given lattice conditions

        Parameters
        ----------
        atom_a, atom_b : Atom objects
            The atoms to be compared
        Returns
        -------
        bonded_bool : bool
            True if the atoms are bonded and False if not
        """
        bonded_bool = atom_a.per_dist(atom_b, self.vectors, ref=self.bonding) <= self.thresh
        return bonded_bool

    # list-y behaviour
    def append(self, element):
        self.atoms.append(element)

    def extend(self, other_mol):
        self.atoms.extend(other_mol.atoms)

    def insert(self, i, element):
        self.atoms.insert(i, element)

    def remove(self, element):
        self.atoms.remove(element)

    def index(self, element):
        return self.atoms.index(element)

    def pop(self, i=-1):
        return self.atoms.pop(i)

    def clear(self):
        self.atoms.clear()

    def count(self, element):
        return self.atoms.count()

    def __add__(self, other_mol):
        try_ismol(other_mol)
        return Mol(deepcopy(self).atoms + other_mol.atoms, vectors = self.vectors, bonding = self.bonding, thresh = self.thresh)

    def __len__(self):
        return len(self.atoms)

    def __eq__(self, other):
        return self.atoms == other.atoms

    def __getitem__(self, index):
        return self.atoms[index]

    def __setitem__(self, index, value):
        self.atoms[index] = value
        return

    def __contains__(self, elem):
        return self.atoms.__contains__(elem)

    def copy(self):
        return deepcopy(self)

    def write_xyz(self, name):
        """Write an xyz file of the Mol"""
        ef.write_xyz(name, self.atoms)

    def empty_mol(self):
        """Return an empty mol with the same properties"""
        new_mol = deepcopy(self)
        new_mol.atoms = []
        return new_mol

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

        selected = self.copy()
        selected.atoms = deepcopy([self[i] for i in labels])
        remaining = self.copy()
        for atom in selected:
            if atom in remaining:
                remaining.remove(atom)

        old_atoms = selected.copy()

        # While there are atoms to add
        cont = True
        while cont:
            cont = False
            new_atoms = Mol([])
            for old in old_atoms:
                tmp_remaining = remaining.copy()
                for rem in remaining:
                    if self.bonded(old,rem):
                        new_atoms.append(rem)
                        selected.append(rem)
                        tmp_remaining.remove(rem)
                        cont = True  # An atom was added so continue loop
                remaining = tmp_remaining
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
        selected_old = self.copy()
        selected_old.atoms = [self[i] for i in labels]

        # Mol of selected atoms where the periodic image
        # atoms are translated back to form a molecule
        selected_img = selected_old.copy()

        remaining = self.copy()
        for atom in selected_old:
            if atom in remaining:
                remaining.remove(atom)

        old_atoms = selected_old.copy()

        # While there are atoms to add
        cont = True
        while cont == True:
            cont = False
            new_atoms = Mol([])
            for old in old_atoms:
                tmp_remaining = remaining.copy()
                for rem in remaining:
                    # contains the distance from the point or image and the
                    # coordinates of the point or image
                    dist, per_img = old.per_dist(rem, self.vectors, ref=self.bonding, new_pos=True)
                    # if the atom is close enough to be part of the molecule
                    if dist <= self.thresh:
                        new_atoms.append(per_img)
                        selected_old.append(rem)
                        selected_img.append(per_img)
                        tmp_remaining.remove(rem)
                        cont = True  # An atom was added so continue loop
                remaining = tmp_remaining
                old_atoms = new_atoms


        if old_pos:
            return selected_img, selected_old
        else:
            return selected_img

    def segregate(self):
        """Separate current Mol in a list of Mols of different molecules"""
        molecules = []  # list of molecules
        remaining = self.copy()

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
        new_cell_atoms = deepcopy(
            [a for a in self.atoms if a not in scattered_mol])
        new_cell = self.copy()
        new_cell.atoms = new_cell_atoms

        for atom in new_mol:
            new_cell.append(atom.copy())
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
        remaining = self.copy()

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

    def supercell(self, trans):
        """
        Return a supercell of I x J x K

        Parameters
        ----------
        trans : array-like of length 3
            Multiplications of the primitive cell
        Returns
        -------
        supercell : Mol object
            New supercell with adjusted lattice vectors

        """
        # make the input into a np array
        trans = np.array(trans)

        new_cell = self.empty_mol()
        for a_mult in range(trans[0]):
            for b_mult in range(trans[1]):
                for c_mult in range(trans[2]):
                    vector = a_mult * \
                        self.vectors[0] + b_mult * \
                        self.vectors[1] + c_mult * self.vectors[2]
                    new_atoms = Mol([i.v_translated(vector)
                                     for i in self.atoms])
                    new_cell += new_atoms
        out_vec = (self.vectors.T * trans.transpose()).T
        new_cell.vectors = out_vec
        return new_cell

    def centered_supercell(self, trans, from_origin=False):
        """
        Make a bigger supercell out of an input cell.

        The cell is multiplied positively and negatively through each lattice
        vector so that the supercluster ends up being
        (1+2*trans[0])*(1+2*trans[1])*(1+2*trans[2]) times larger. For example if the
        input is 1,1,1 for a cubic unit cell, the output will be the original unit
        cell surrounded by 26 other unit cells forming a total 3x3x3 cube.

        Alternatively, the multiplication can be centered around the origin, a corner of the
        unit cell, instead of the centre. In that case the supercluster ends up being
        only (2*trans[0])*(2*trans[1])*(2*trans[2])

        Parameters
        ----------
        trans : numpy array of length 3
            Multiplications of the primitive cell
        from_origin : bool
            Determines the kind of multiplication. True is corner of the cell as
            the center, False is middle of the cell.

        Returns
        -------
        mega_cell : Mol object
            The resulting supercell

        """
        trans_series = [0, 0, 0]
        for i, tra in enumerate(trans):
            if from_origin:
                trans_series[i] = list(range(-tra, tra))
            else:
                trans_series[i] = list(range(-tra, tra + 1))
        trans_series = np.array(trans_series)

        new_cell = self.empty_mol()
        for a_mult in trans_series[0]:
            for b_mult in trans_series[1]:
                for c_mult in trans_series[2]:
                    vector = a_mult * \
                        self.vectors[0] + b_mult * \
                        self.vectors[1] + c_mult * self.vectors[2]
                    new_atoms = Mol([i.v_translated(vector)
                                     for i in self.atoms])
                    new_cell += new_atoms
        out_vec = (self.vectors.T * trans.transpose()).T
        new_cell.vectors = out_vec
        return new_cell

    def trans_from_rad(self, clust_rad):
        """
        Generate the translations necessary to encapsulate a sphere of given rad

        Parameters
        ----------
        clust_rad : float
            Radius defining a sphere

        Returns
        -------
        trans_count : 3 x 1 numpy array
            The translations required for the unit cell to contain the sphere

        """

        # determine how many unit cells we need
        vectors = deepcopy(self.vectors)

        # vectors normal to faces
        a_perp = np.cross(vectors[1], vectors[2])
        b_perp = np.cross(vectors[2], vectors[0])
        c_perp = np.cross(vectors[0], vectors[1])

        # the three normalised unit vectors
        perp = np.array([a_perp / np.linalg.norm(a_perp), b_perp /
                         np.linalg.norm(b_perp), c_perp / np.linalg.norm(c_perp)])

        trans_count = np.array([1, 1, 1])

        # distances from faces
        distances = np.array([0.0, 0.0, 0.0])

        new_vectors = deepcopy(vectors)

        for comp in range(3):
            while True:
                trans_count[comp] += 1
                distances[comp] = np.dot(new_vectors[comp], perp[comp])
                new_vectors[comp] = trans_count[comp] * vectors[comp]
                if distances[comp] > clust_rad:
                    break
        trans_count -= np.array([1, 1, 1])
        return trans_count

    def make_cluster(self, clust_rad, mode = 'exc'):
        """
        Generate a cluster of molecules from a primitive cell

        This first makes a supercell of the correct size which will contain with
        one additional buffer shell. Then the sphere is generated from this new
        supercell by connectivity.

        Parameters
        ----------
        clust_rad : float
            Radius defining a sphere. All molecules with atoms in the sphere are
            to be grabbed
        inc : mode
            Switches between inclusive and exclusive selecting. Inclusive,
            'inc', selects all molecules which have atoms within the radius.
            Exclusive, 'exc', selects all molecules fully in the radius.
            Default: false
        Returns
        -------
        cluster : Mol object
            Spherical cluster of molecules from their crystal positions

        """
        trans = self.trans_from_rad(clust_rad)
        if mode == 'inc':
            trans += np.array([1,1,1]) # one buffer cell layer
        supercell = self.centered_supercell(trans, from_origin=True)
        # atoms within the sphere of rad clust_rad
        seed_atoms = Mol([])
        for atom in supercell:
            if atom.v_dist([0, 0, 0]) < clust_rad:
                seed_atoms.append(atom)
        max_mol_len = 0
        if mode == 'exc':
            while len(seed_atoms) > 0:
                mol = seed_atoms.select(0)
                if len(mol) > max_mol_len:
                    max_mol_len = len(mol)
                    clust_atoms = Mol([])
                if len(mol) == max_mol_len:
                    clust_atoms += mol
                for atom in mol:
                    seed_atoms.remove(atom)
        if mode == 'inc':
            clust_atoms = Mol([])
            max_mol_len = len(supercell.select(supercell.index(seed_atoms[0])))

            while len(seed_atoms) > 0:
                print(len(seed_atoms))
                mol_tmp = seed_atoms.select(0) # The part of the mol detected in seed_atoms
                if len(mol_tmp) < max_mol_len:
                    # The whole mol, which could potentially include even more seed_atoms
                    mol = supercell.select(supercell.index(seed_atoms[0]))
                else:
                    mol = mol_tmp
                clust_atoms += mol
                for atom in mol_tmp:
                    seed_atoms.remove(atom)
                for atom in mol:
                    supercell.remove(atom)
                    # remove all atoms of the mol which are part of seed_atoms
                    try:
                        seed_atoms.remove(atom)
                    except ValueError:
                        pass

        return clust_atoms

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
        """Move all atoms to fractional coordinates"""

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

    def confined(self):
        """Move all atoms to fit inside the primitive cell"""
        frac_mol = self.dir_to_frac_pos()
        out_mol = frac_mol.frac_to_dir_pos()

        return out_mol

    def centered_mols(self, labels, return_trans = False):
        """
        Return the molecules translated at the origin with a corresponding cell

        Parameters
        ----------
        labels : int or list of ints
            The labels of the atoms to select
        print_centro : bool
            Print the translation vector which was detected as -centroid
        Returns
        -------
        mol : Mol object
            The selected molecules with their centroid at the origin
        mod_cell : Mol object
            The new confined cell corresponding to the now translated molecules

        """
        mol, mod_cell = self.complete_mol(labels)
        centro = mol.centroid()
        mol.translate(-centro)
        mod_cell.translate(-centro)
        mod_cell = mod_cell.confined()

        if return_trans:
            return mol, mod_cell, -centro
        else:
            return mol, mod_cell


    def es_pot(self, position):
        """
        Return the electorstatic potential generated by this Mol

        Parameters
        ----------
        position : 3x1 np array
            The point at which the potential should be evaluated
        Returns
        -------
        tot_pot : float
            The total potential

        """
        tot_pot = 0
        for atom in self:
            tot_pot += atom.es_pot(position)
        return tot_pot

    def change_charges(self, charges):
        """
        Change all of the charges of the constituent atoms at once

        Parameters
        ----------
        charges : array-like of floats
            Contains all of the new charges. IMPORTANT: they need to be in the
            order corresponding to self.atoms

        """
        for i, atom in enumerate(self.atoms):
            atom.q = charges[i]
        return

    def charges(self):
        """Return an array of charges"""
        l_char = []
        for atom in self.atoms:
            l_char.append(atom.q)
        arr_char = np.array(l_char)
        return arr_char

    def raw_assign_charges(self, charges):
        """Assign the charges from an array-like to the atoms"""
        for char,at in zip(charges,self.atoms):
            at.q=char
        return

    def populate(self, reference_mol):
        """
        Assign charges to the Mol by comparing to the connectivity of a
        reference

        Parameters
        ----------
        reference_mol : Mol object
            Charged molecule or cell

        """
        # This is a naughty in-function import to prevent a circular dependency.
        # The reason is that assign_charges functions are grouped up with the
        # executable script which needs to read_file and in turn use mol.py
        # Some careful refactoring should fix this
        import fromage.scripts.assign_charges as ac
        ac.assign_charges(reference_mol, self)
        pass
