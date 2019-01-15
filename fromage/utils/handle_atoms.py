"""Some functions to handle lists of atoms
"""
import numpy as np
from copy import copy

from fromage.utils.atom import Atom


def select(max_r, atoms, label):
    """
    Select a molecule out of a list of atoms.

    The function puts a labelled atom in a list of selected atoms. It then
    checks for nearby atoms to the selected atoms. If any are found and are not
    already selected, they are added to the list. The loop ends when it is
    executed and not atoms are added.

    Parameters
    ----------
    max_r : float
        Maximum distance which is considered a bond
    atoms : list of Atom objects
        Cluster of atoms which from which a molecule needs to be extracted
    label : int
        The number of the atom from which the molecule is generated. Note that
        it is the Python label, so starting from 0, unlike the label in most
        chemistry software which counts from 1

    Returns
    -------
    selected : list of Atom objects
        The atoms belonging to the molecule which is selected

    """

    new_atoms = [atoms[label]]
    selected = new_atoms
    n = True
    while n == True:
        n = False
        old_atoms = new_atoms
        new_atoms = []
        for i in old_atoms:
            for j in atoms:
                # find atoms bonded to the newest atoms in the list
                if i.c_dist(j.x, j.y, j.z) <= max_r and j not in selected:
                    new_atoms.append(j)
                    selected.append(j)
                    # if no more atoms are found stop the loop
                    n = True

    return selected


def select_per(max_r, atoms, label, vectors):
    """
    Select a molecule out of a list of atoms in a periodic system.

    Parameters
    ----------
    max_r : float
        Maximum distance which is considered a bond
    atoms : list of Atom objects
        Cluster of atoms which from which a molecule needs to be extracted
    label : int
        The number of the atom from which the molecule is generated. Note that
        it is the Python label, so starting from 0, unlike the label in most
        chemistry software which counts from 1
    vectors : 3x3 matrix
        Lattice vectors

    Returns
    -------
    selected : list of Atom objects
        The atoms belonging to the molecule which is selected
    selected_img : list of Atom objects
        The atoms belonging to the molecule which is selected with certain
        atoms translated so that the molecule is fully connected without
        periodic boundaries

    """
    # list of selected atoms from the unit cell
    selected = [atoms[label]]
    # list of selected atoms where the periodic image
    # atoms are translated back to form a molecule
    selected_img = [atoms[label]]

    # n is true as long as there are still molecules to add to the list
    n = True
    while n == True:
        n = False
        for i in selected_img:
            for j in atoms:
                # contains the distance from the point or image and the
                # coordinates of the point or image
                gamma = i.dist_lat(j.x, j.y, j.z, vectors[
                    0], vectors[1], vectors[2])

                # if the atom is close enough to be part of the molecule
                # and is not already part of the molecule
                if gamma[0] <= max_r and j not in selected:
                    selected.append(j)
                    k = copy(j)
                    k.x, k.y, k.z = gamma[1:]
                    selected_img.append(k)
                    n = True

    return selected, selected_img


def multi_select(max_r, atoms, labels, vectors):
    """
    Select multiple molecules in a periodic system

    Parameters
    ----------
    max_r : float
        Maximum distance which is considered a bond
    atoms : list of Atom objects
        Cluster of atoms which from which molecules needs to be extracted
    labels : list of ints
        The numbers of the atom from which the molecules are generated. Note
        that it is the Python label, so starting from 0, unlike the label in
        most chemistry software which counts from 1
    vectors : 3x3 matrix
        Lattice vectors

    Returns
    -------
    selected : list of Atom objects
        The atoms belonging to the molecules which are selected
    selected_img : list of Atom objects
        The atoms belonging to the molecules which are selected with certain
        atoms translated so that the molecules are fully connected without
        periodic boundaries

    """
    selected_mols = []
    selected_img_mols = []

    for label in labels:
        selected_mol, selected_img_mol = select_per(
            max_r, atoms, label, vectors)
        selected_mols.append(selected_mol)
        selected_img_mols.append(selected_img_mol)

        # checks to see if the user selected the same molecule twice
    for moleculeI in selected_mols:
        for moleculeJ in selected_mols:
            if moleculeI != moleculeJ:
                if moleculeI[0] in moleculeJ:
                    raise ValueError(
                        "You have selected several atoms in the same molecule!")
    tmp_a = [item for sublist in selected_mols for item in sublist]
    tmp_b = [item for sublist in selected_img_mols for item in sublist]

    selected_mols = tmp_a
    selected_img_mols = tmp_b

    return selected_mols, selected_img_mols


def make_molecules(atoms, bl):
    """
    Generate list of molecules based on bond length bl

    Parameters
    ----------
    atoms: list of Atom objects
    bl: float
        Bond length in unit of input file
    Returns
    -------
    molecules: list of lists
        List L of length M molecules, where each member of L is a list of atom objects
    """
    molecules = []  # list of molecules
    max_length = 0  # number of atoms in a molecule
    for i, atom in enumerate(atoms):
        if atom not in [val for sublist in molecules for val in sublist]:
            # creates a molecule if atom not already in molecules list
            molecule = select(bl, atoms[i:], 0)
            if len(molecule) > max_length:
                max_length = len(molecule)
                molecules = []
                molecules.append(molecule)
            elif len(molecule) == max_length:
                molecules.append(molecule)
    return molecules


def complete_mol(max_r, atoms, label, vectors):
    """
    Take a cell and complete one molecule.

    The objective is to end up with a unit cell where the molecule of interest
    is complete. The rest of the atoms of the cell can remain intact. Note that
    the input atoms are transformed and are the same as are present in the
    output.

    Parameters
    ----------
    max_r : float
        Maximum distance which is considered a bond
    atoms : list of Atom objects
        Cluster of atoms which from which a molecule needs to be extracted
    label : list of ints
        The number of the atom from which the molecule is generated. Note that
        it is the Python label, so starting from 0, unlike the label in most
        chemistry software which counts from 1
    vectors : 3x3 matrix
        Lattice vectors

    Returns
    -------
    full_mol_trans : list of Atom objects
        The now complete molecule
    atoms : list of Atom objects
        The cell with the completed molecule. NB the Atom objects of the completed molecule
        reference the same objects as full_mol_trans

    """
    # the atoms which are part of the same selected molecule.
    # first with intact coordinates and second with appropriate
    # atoms translated to join up the molecule with their image
    full_mol, full_mol_trans = multi_select(max_r, atoms, label, vectors)

    # using an explicit loop does not work for some reason,
    # use list comprehension
    atoms = [a for a in atoms if a not in full_mol]

    for atom in full_mol_trans:
        atoms.append(atom)
    return full_mol_trans, atoms


def complete_cell(atoms, vectors, max_bl=1.7):
    """
    Return a cell where atoms have been translated to complete all molecules of
    the cell

    Parameters
    ----------
    atoms : list of Atom objects
        The input truncated unit cell
    vectors : 3 x 3 numpy array
        The unit cell vectors
    max_bl : float
        The maximum interatomic distance which is considered a bond
    Returns
    -------
    out_cell : list of Atom objects
        The new untruncated cell
    full_mol_l : list of lists of Atom objects
        Each molecule in the untruncated cell

    """
    full_mol_l = []

    while len(atoms) != 0:
        full_mol, cell = complete_mol(max_bl, atoms, [0], vectors)
        full_mol_l.append(full_mol)
        atoms = cell
        for atom in full_mol:
            if atom in atoms:
                atoms.remove(atom)

    out_cell = []
    for mol in full_mol_l:
        out_cell.extend(mol)
    return out_cell, full_mol_l


def find_centroid(atoms):
    """Find the centroid of a list of atoms."""
    N = len(atoms)
    baryX, baryY, baryZ = 0, 0, 0
    for atom in atoms:
        baryX += atom.x / N
        baryY += atom.y / N
        baryZ += atom.z / N
    return (baryX, baryY, baryZ)


def center_mol(atoms):
    """Translate molecules to center"""
    c_x, c_y, c_z = find_centroid(atoms)
    for atom in atoms:
        atom.translate(-c_x, -c_y, -c_z)

    return


def translate_atoms(atoms, vector):
    """
    Translate some atoms by a vector

    Parameters
    ----------
    atoms : list of Atom objects
        Untranslated atoms
    vector : 3 x 1 numpy array
        Translation vector

    """
    for atom in atoms:
        atom.v_translate(vector)
    return


def supercell(atoms, vectors, trans):
    """
    Make a supercell of I x J x K

    Parameters
    ----------
    atoms : list of Atom objects
        Input cell atoms
    vectors : 3 x 3 numpy array
        Lattice vectors
    trans : numpy array of length 3
        Multiplications of the primitive cell
    Returns
    -------
    supercell : list of Atom objects
        New supercell of Atoms
    out_vec : 3 x 3 numpy array
        The lattice vectors of the new cell

    """
    cart = [0, 1, 2]
    supercell = copy(atoms)
    for comp in cart:
        if trans[comp] != 1:
            for mult in range(trans[comp])[1:]:
                new_atoms = [i.v_translated(mult * vectors[comp])
                             for i in supercell]
                supercell += new_atoms
    out_vec = (vectors.T * trans.transpose()).T
    return supercell, out_vec


def make_mega_cell(atoms, traAN, traBN, traCN, vectors, alt_multi=False):
    """
    Make a bigger supercell out of an input cell.

    The cell is multiplied positively and negatively through each lattice
    vector so that the supercluster ends up being
    (1+2*traAN)*(1+2*traBN)*(1+2*traCN) times larger. For example if the
    input is 1,1,1 for a cubic unit cell, the output will be the original unit
    cell surrounded by 26 other unit cells forming a total 3x3x3 cube.

    Alternatively, the multiplication can be centered around a corner of the
    unit cell instead of the centre. In that case the supercluster ends up being
    only (2*traAN)*(2*traBN)*(2*traCN)

    Parameters
    ----------
    atoms : list of Atom objects
        Atoms in the cell
    traAN,traBN,traCN : ints
        Number of times the cell should be multiplied in each cell direction,
        positively and negatively
    vectors : 3x3 matrix
        Lattice vectors
    alt_multi : bool
        Determines the kind of multiplication. True is centre of the cell as
        the origin, False is corner of the cell.
    Returns
    -------
    mega_cell : list of Atom objects
        The resulting supercell

    """
    mega_cell = []

    if alt_multi:
        traA = range(-traAN, traAN)
        traB = range(-traBN, traBN)
        traC = range(-traCN, traCN)
    else:
        traA = range(-traAN, traAN + 1)  # +1 because that is how range() works
        traB = range(-traBN, traBN + 1)
        traC = range(-traCN, traCN + 1)

    # multiplying the cell
    for i in traA:
        for j in traB:
            for k in traC:
                traX = vectors[0][0] * i + \
                    vectors[1][0] * j + vectors[2][0] * k
                traY = vectors[0][1] * i + \
                    vectors[1][1] * j + vectors[2][1] * k
                traZ = vectors[0][2] * i + \
                    vectors[1][2] * j + vectors[2][2] * k

                for atom in atoms:
                    mega_cell.append(atom.translated(traX, traY, traZ))
    return mega_cell


def make_cluster(atoms, clust_rad, max_bl):
    """
    Generate a cluster of molecules from a cluster of atoms.

    Make sure the input cluster of atoms is large enough to contain the
    desired cluster of molecules. You can generate a large cluster of atoms with
    make_mega_cell. The molecules are defined as any with atoms falling within
    a certain radius of the origin.

    Parameters
    ----------
    atoms : list of Atom objects
        A large cluster of atoms
    clust_rad : float
        Radius determining the initial atoms of the cluster of molecules. These
        atoms are then the seed for full molecules
    max_bl : float
        Maximum distance which counts as a bond

    Returns
    -------
    clust_atoms : list of Atom objects
        A list of atoms which form a cluster of molecules

    """
    # atoms within the sphere of rad clust_rad
    seed_atoms = []

    for atom in atoms:
        if atom.c_dist(0, 0, 0) < clust_rad:
            seed_atoms.append(atom)

    # atoms in the cluster (seed_atoms + atoms to complete molecules)
    clust_atoms = []
    for i, atom in enumerate(seed_atoms):
        if atom not in clust_atoms:
            mol2Add = select(max_bl, atoms, atoms.index(atom))
            for atom2Add in mol2Add:
                clust_atoms.append(atom2Add)
    return clust_atoms


def array2atom(template, pos):
    """
    Turn an array of the form x1, y1, z1, x2, y2, z2 etc. into a list of Atom
    objects

    Parameters
    ----------
    template : list of Atom objects
        A list of the same length of the desired one used to determine the
        elements of the atoms
    pos : list of floats
        List of coordinates of the form x1, y1, z1, x2, y2, z2 etc.
    Returns
    -------
    out_atoms : list of Atom objects
        Resulting atoms

    """
    sliced_pos = [pos[i:i + 3] for i in range(0, len(pos), 3)]
    out_atoms = []
    for atom in zip(template, sliced_pos):
        new_atom = Atom(atom[0].elem, atom[1][0], atom[1][1], atom[1][2], 0)
        out_atoms.append(new_atom)
    return out_atoms
