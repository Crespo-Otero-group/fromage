default_thresh = {'dis': 1.8,
                  'cov': 0.2,
                  'vdw': -0.3}


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
        raise TypeError("Unrecognised bonding type: " + bonding)
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
    thresh_str = ''
    # check if bonding has been specified
    for i_bonding in bondings:
        if i_bonding in in_str:
            bonding = i_bonding
    # if there is bonding, try to find threshold
    if bonding:
        stripped = in_str.replace(bonding, '')
        # if there is still a thresh in this string
        if stripped:
            thresh_str = stripped
    # if only the thresh is specified
    elif in_str:
        thresh_str = in_str
    # if both present
    if bonding and thresh_str:
        self.set_bonding(bonding=bonding, thresh=float(thresh_str))
    # if only bonding
    if bonding and not thresh_str:
        self.set_bonding(bonding=bonding)
    # if only thresh
    if thresh_str and not bonding:
        self.set_bonding(thresh=float(thresh_str))
    if not thresh_str and not bonding:
        self.set_bonding()
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
    bonded_bool = atom_a.per_dist(
        atom_b, self.vectors, ref=self.bonding) <= self.thresh
    return bonded_bool
