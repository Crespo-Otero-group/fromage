from copy import deepcopy

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
    import fromage.utils.mol as mol_init

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
        new_atoms = mol_init.Mol([])
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
    import fromage.utils.mol as mol_init

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
        new_atoms = mol_init.Mol([])
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
