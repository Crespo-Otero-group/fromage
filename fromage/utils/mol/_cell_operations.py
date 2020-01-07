from copy import deepcopy
import numpy as np


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
    import fromage.utils.mol as mol_init
    # make the input into a np array
    trans = np.array(trans)

    new_cell = self.empty_mol()
    for a_mult in range(trans[0]):
        for b_mult in range(trans[1]):
            for c_mult in range(trans[2]):
                vector = a_mult * \
                    self.vectors[0] + b_mult * \
                    self.vectors[1] + c_mult * self.vectors[2]
                new_atoms = mol_init.Mol([i.v_translated(vector)
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
    import fromage.utils.mol as mol_init

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
                new_atoms = mol_init.Mol([i.v_translated(vector)
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

def supercell_for_cluster(self, clust_rad, mode='exc', central_mol=None):
    """
    Make a supercell which will be used to make a cluster

    Parameters
    ----------
    clust_rad : float
        The radius of the cluster
    mode : str
        'exc' is for exclusive clusters, whereas 'inc' is for inclusive clusters
        which will need an extra layer of unit cells
    central_mol : Mol
        Molecule which serves as a center for the cluster (optional)
    Returns
    -------
    out_supercell : Mol
        The supercell from which the cluster will be taken

    """
    # if there is a central mol, account for nearest neighbour molecules
    # bleeding out of the original radius
    if central_mol:
        central_rad = 0
        for atom in central_mol:
            dis = atom.v_dist([0, 0, 0])
            if dis < central_rad:
                central_rad = dis
    # get the translations of the unit cell necessary to enclose the required mols
        trans = self.trans_from_rad(clust_rad + central_rad)
    else:
        trans = self.trans_from_rad(clust_rad)

    # if the cluster is inclusive, then extra mols might be required from
    # an additional layer of the supercell
    if mode == 'inc':
        trans += np.array([1, 1, 1])  # one buffer cell layer

    # make a supercell which includes the desired cluster
    out_supercell = self.centered_supercell(trans, from_origin=True)

    return out_supercell

def gen_exclusive_clust(self, seed_atoms):
    """
    Remove all non complete molecules

    This only works if the input contains at least one full molecule

    Parameters
    ----------
    seed_atoms : Mol
        Aggregate of molecules, not necessarily all complete
    Returns
    -------
    out_clust : Mol
        The full molecules of seed_atoms

    """
    import fromage.utils.mol as mol_init

    max_mol_len = 0
    while len(seed_atoms) > 0:
        # pick out a molecule from the seed atoms
        mol = seed_atoms.select(0)

        # if the molecule is the biggest so far
        if len(mol) > max_mol_len:
            # molecules are supposed to be this long now
            max_mol_len = len(mol)
            out_clust = mol_init.Mol([])
        # if the molecule is the right size
        if len(mol) == max_mol_len:
            # include this molecule
            out_clust += mol

        # discard the molecule from seed atoms
        for atom in mol:
            seed_atoms.remove(atom)

    return out_clust


def gen_inclusive_clust(self, seed_atoms, supercell):
    """
    Select all complete molecules of supercell which contain seed atoms

    This only works if the input contains at least one full molecule

    Parameters
    ----------
    seed_atoms : Mol
        Aggregate of molecules, not necessarily all complete
    supercell : Mol
        Supercell which contains all seed atoms. It should have at least one
        buffer layer of unit cells around the seed atoms
    Returns
    -------
    out_clust : Mol
        The full molecules of supercell which contain seed atoms

    """
    import fromage.utils.mol as mol_init

    max_mol_len = 0

    out_clust = mol_init.Mol([])
    # here, the molecule with the atom seed_atoms[0] is necessarily complete
    # in supercell
    max_mol_len = len(supercell.select(supercell.index(seed_atoms[0])))

    while len(seed_atoms) > 0:
        # the part of the mol detected in seed_atoms
        mol_tmp = seed_atoms.select(0)
        if len(mol_tmp) < max_mol_len:
            # The whole mol, which could potentially include even more
            # seed_atoms
            mol = supercell.select(supercell.index(seed_atoms[0]))
        else:
            mol = mol_tmp
        out_clust += mol
        for atom in mol_tmp:
            seed_atoms.remove(atom)
        for atom in mol:
            supercell.remove(atom)
            # remove all atoms of the mol which are part of seed_atoms
            try:
                seed_atoms.remove(atom)
            except ValueError:
                pass

    return out_clust

def make_cluster(self, clust_rad, mode='exc', central_mol=None):
    """
    Generate a cluster of molecules from a primitive cell

    This first makes a supercell of the correct size which will contain with
    one additional buffer shell. Then the sphere is generated from this new
    supercell by connectivity.

    A central molecule can also be supplied which will turn the spheres
    defining the clusters into the union of spheres stemming from each atom
    of the central molecule.

    This algorithm is designed for crystals where one same molecule does not
    extend into more than two unit cells in the case of inclusive clusters.

    Parameters
    ----------
    clust_rad : float
        Radius for generating the cluster. If no central molecule is specified,
        this will generate seed atoms in a sphere from the radius
    mode : str
        Switches between inclusive and exclusive selecting. Inclusive,
        'inc', selects all molecules which have atoms within the radius.
        Exclusive, 'exc', selects all molecules fully in the radius.
        Default: 'exc'
    central_mol : Mol
        If this is supplied, the central molecule will act as a kernel for the
        cluster which will end up being of the appropriate shape. (optional)
    Returns
    -------
    cluster : Mol object
        Cluster of molecules from their crystal positions

    """
    import fromage.utils.mol as mol_init
    # generate a supercell which will include the cluster.
    # inclusive clusters will have an extra layer of supercell.
    # if a central mol is supplied, the supercell will include the whole
    # molecule and the supplied radius.
    supercell = self.supercell_for_cluster(clust_rad, mode=mode, central_mol=central_mol)

    # seed_atoms will initialise the cluster
    seed_atoms = mol_init.Mol([])

    # get seed atoms in the shape of the central mol if pertinent
    if central_mol:
        for atom_i in supercell:
            for atom_j in central_mol:
                if atom_i.dist(atom_j) < clust_rad:
                    seed_atoms.append(atom_i)
                    break

    # get spherical seedatoms otherwise
    else:
        for atom in supercell:
            if atom.v_dist([0, 0, 0]) < clust_rad:
                seed_atoms.append(atom)

    # refine seed atoms
    max_mol_len = 0
    # remove incomplete molecules
    if mode == 'exc':
        clust_atoms = self.gen_exclusive_clust(seed_atoms)
    # complete incomplete molecules
    elif mode == 'inc':
        clust_atoms = self.gen_inclusive_clust(seed_atoms, supercell)
    else:
        raise ValueError("Invalid cluster generation mode. Use 'exc' or 'inc'")

    return clust_atoms


def centered_mols(self, labels, return_trans=False):
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

def confined(self):
    """Move all atoms to fit inside the primitive cell"""
    frac_mol = self.dir_to_frac_pos()
    out_mol = frac_mol.frac_to_dir_pos()

    return out_mol
