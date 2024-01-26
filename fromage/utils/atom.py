"""Defines the Atom object"""

import numpy as np
from collections import Counter
from copy import deepcopy
from fromage.utils import per_table as per
from fromage.fdist import _fdist as fd


class Atom(object):
    """
    Object representing an atom.

    Sometimes also used to represent point charges as atoms of element "point".
    Several functions are present like translate or find_centroid.

    Attributes
    ----------
    x,y,z : floats
        Cartesian coordinates
    q : float
        Partial atomic charge
    connectivity : frozenset of tuples
        The set is ((atom kind,connectivity order),amount) and is set via a
        function which takes the connectivity matrix as argument
    kind : tuple
        Tuple of (atom element,connectivity). This defines the kind of atom
    total_e : int
        Atomic number
    valence : int
        Number of valence electrons
    vdw : float
        Van der Waals radius in Angstrom
    cov : float
        Covalent radius in Angstrom
    """

    def __init__(self, elemIn="H", xIn=0.0, yIn=0.0, zIn=0.0, qIn=0.0, num=1):
        self.elem = elemIn
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.q = 0.0
        self.num = 1
        # Atom objects with no charge can have feel a finite electostatic
        # potential which we include as
        self.es = 0.0
        # The connectivity is a frozenset (because a list would have a built-in ordering)
        # of the tuples of the form (A,N) where A is an tuple of different
        # distances to atoms e.g. ("C",4) if there is a carbon 4 bonds away. N
        # is the amount of carbons 4 bonds away
        self.connectivity = None
        # Kind is a tuple of (elem,connectivity) and as such is enough to define
        # an atom type at least as well as it would be defined in a forcefield
        # e.g. in acrolein: This is a C atom with 1 O 1-away, an H 1-away, a C
        # 1-away, an H 2-away, a C 2-away and 2 H 3-away
        self.kind = None

        # deal with some sneaky int that may be disguised as float
        try:
            self.x = float(xIn)
            self.y = float(yIn)
            self.z = float(zIn)
            self.q = float(qIn)

        except ValueError:
            print("Some coordinates or charges cannot be cast to float!")

        table = per.periodic
        self.at_num = table[self.elem.lower()]["at_num"]
        self.valence_e = table[self.elem.lower()]["valence_e"]
        self.cov = table[self.elem.lower()]["cov"]
        self.vdw = table[self.elem.lower()]["vdw"]
        self.mass = table[self.elem.lower()]["mass"]
        # to string methods to be used mainly for debugging and .qc file

    def __repr__(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(self.elem, self.x, self.y, self.z, self.q)

    def __str__(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(self.elem, self.x, self.y, self.z, self.q)

        # equality function
    def __eq__(self, other):
        return self.elem.lower() == other.elem.lower() and self.dist(other) < 1e-5 and self.q - other.q < 1e-5
#        return self.elem.lower() == other.elem.lower() and self.x == other.x and self.y == other.y and self.z == other.z and self.q == other.q

    def copy(self):
        return deepcopy(self)
    def set_pos(self, pos_array):
        """
        Assign coordinates via numpy array

        Parameters
        ----------
        pos_array : 3 x 1 list-like
            The position to assign to the atom

        """

        self.x = pos_array[0]
        self.y = pos_array[1]
        self.z = pos_array[2]

        return

    def get_pos(self):
        """Return np array of coord"""
        out_arr = np.array([self.x, self.y, self.z])

        return out_arr

    def very_close(self, other, thresh=0.001):
        """Check if two atoms are very close together"""
        x_cond = abs(self.x - other.x) < thresh
        y_cond = abs(self.y - other.y) < thresh
        z_cond = abs(self.z - other.z) < thresh

        cond = x_cond and y_cond and z_cond
        return cond

    def xyz_str(self):
        """Return a string of the atom in xyz format"""
        return "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(self.elem, self.x, self.y, self.z)

    def c_dist2(self, x1, y1, z1):
        """Return distance squared of the atom from a point"""
        # Use for no C++ version
#        r = (self.x - x1) ** 2 + (self.y - y1) ** 2 + (self.z - z1) ** 2
        r = fd.dist2(self.x, self.y, self.z, x1, y1, z1)
        return r

    def c_dist(self, x1, y1, z1):
        """Return distance of the atom from a point"""
        #r = np.sqrt(self.dist2(x1, y1, z1))
        r = fd.dist(self.x, self.y, self.z, x1, y1, z1)
        return r

    def v_dist(self, position):
        """Return distance of the atom from a point defined by an array-like"""
        r = self.c_dist(position[0], position[1], position[2])
        return r

    def v_dist2(self, position):
        """Return distance of the atom from a point defined by an array-like"""
        r = self.c_dist2(position[0], position[1], position[2])
        return r

    def dist(self, other_atom, ref='dis'):
        """Return a distance between atoms starting form different ref points"""
        distance_types = {'dis' : self.dist_dis,
                        'cov' : self.dist_cov,
                        'vdw' : self.dist_vdw}

        r = distance_types[ref](other_atom)
        return r

    def dist_dis(self, other_atom):
        """Return the distance between centres of atoms"""
        r = self.c_dist(other_atom.x,other_atom.y,other_atom.z)
        return r

    def dist_cov(self, other_atom):
        """Return the distance between covalent spheres of atoms"""
        r = self.dist_dis(other_atom) - self.cov - other_atom.cov
        return r

    def dist_vdw(self, other_atom):
        """Return the distance between vdw spheres of atoms"""
        r = self.dist_dis(other_atom) - self.vdw - other_atom.vdw
        return r

    def at_lap(self, other_atom):
        """Return the overlap between vdw radii"""
        r = self.vdw + other_atom.vdw - self.dist_dis(other_atom)
        return r

    def dist_lat(self, x1, y1, z1, aVec, bVec, cVec, order=1, ref='dis'):
        """
        Find the shortest distance to a point in a periodic system.

        Parameters
        ----------
        x1,y1,z1 : floats
            Cartesian coordinates of the target point
        aVec,bVec,cVec : 3x1 array-likes
            Unit cell vectors
        order : positive int
            The amount of translations to be considered. Order 1 considers a
            translation by -1, 0 and 1 of each lattice vector and all resulting
            combination. Order 2 is [-2, -1, 0, 1, 2] and so onzx

        Returns
        -------
        rMin : float
            Minimal distance to the point
        x3,y3,z3 : floats
            Coordinates of the closest image to the point

        """

        in_pos = np.array([x1, y1, z1])
        vectors = np.array([aVec, bVec, cVec])
        multipliers = np.arange(-order, order + 1)

        # sets comprised of the ranges of lattice vector values
        aSet = [i * vectors[0] for i in multipliers]
        bSet = [i * vectors[1] for i in multipliers]
        cSet = [i * vectors[2] for i in multipliers]

        # minimum r distance
        rMin = float("inf")

        # loop over all possible translations of the input point
        for trans1 in aSet:
            for trans2 in bSet:
                for trans3 in cSet:
                    img_pos = in_pos + trans1 + trans2 + trans3
                    # r=np.linalg.norm(self.my_pos-img_pos)
                    r = self.c_dist(img_pos[0], img_pos[1], img_pos[2])
                    # if this particular translation of the point is the closest
                    # to the atom so far
                    if r < rMin:
                        rMin = r
                        # image coordinates
                        x3 = img_pos[0]
                        y3 = img_pos[1]
                        z3 = img_pos[2]
        return rMin, x3, y3, z3

    def per_dist(self, other_atom, vectors, ref='dis', new_pos=False):
        """
        Find the shortest distance to another atom in a periodic system.

        Parameters
        ----------
        other_atom : Atom object
            The atom which to which the distance is being calculated
        vectors : 3 x 3 numpy array
            Unit cell vectors
        order : positive int
            The amount of translations to be considered. Order 1 considers a
            translation by -1, 0 and 1 of each lattice vector and all resulting
            combination. Order 2 is [-2, -1, 0, 1, 2] and so onzx

        Returns
        -------
        r_min : float
            Minimal distance to the point
        at_img : floats (optional)
             Closest image of the atom being targeted

        """
        # First put the other atom inside the cell
        other_atom = other_atom.put_in_cell(vectors)

        multipliers = np.array([-1, 0, 1])

        # sets comprised of the ranges of lattice vector values
        a_set = [i * vectors[0] for i in multipliers]
        b_set = [i * vectors[1] for i in multipliers]
        c_set = [i * vectors[2] for i in multipliers]

        # minimum r distance
        r_min = float("inf")

        # is the minimal distance unique?
        unique = True

        # loop over all possible translations of the input point
        for trans_a in a_set:
            for trans_b in b_set:
                for trans_c in c_set:
                    cell_origin = trans_a + trans_b + trans_c
                    tmp_img_atom = other_atom.v_translated(cell_origin)
                    r = self.dist(tmp_img_atom, ref=ref)
                    if r <= r_min:
                        if r < r_min:
                            unique = True
                        if r == r_min:
                            unique = False
                        r_min = r
                        at_img = tmp_img_atom
        if not unique:
            print("WARNING: the closest periodic image is ill-defined")
        if new_pos:
            return r_min, at_img
        else:
            return r_min

    def put_in_cell(self, vectors):
        """
        Return a new atom at a position inside the parallelepiped cell
        """
        # transpose to get the transformation matrix
        M = np.transpose(vectors)
        # inverse transformation matrix
        U = np.linalg.inv(M)

        dir_pos = np.array([self.x, self.y, self.z])
        frac_pos = np.dot(U, dir_pos)
        # translate to the range [0.1]
        new_frac_pos = [i % 1 for i in frac_pos]

        new_dir_pos = np.dot(M, new_frac_pos)

        new_at = deepcopy(self)
        new_at.x = new_dir_pos[0]
        new_at.y = new_dir_pos[1]
        new_at.z = new_dir_pos[2]

        return new_at


    def translated(self, x1, y1, z1):
        """Return a new atom which is a translated copy."""
        xout, yout, zout = self.x, self.y, self.z
        xout += x1
        yout += y1
        zout += z1
        outAtom = Atom(self.elem, xout, yout, zout, self.q)
        return outAtom

    def v_translated(self, vec_trans):
        """Return a new atom which is a translated copy."""
        old_pos = np.array([self.x, self.y, self.z])
        new_pos = old_pos + vec_trans
        out_atom = Atom(self.elem, new_pos[0], new_pos[1], new_pos[2], self.q)
        out_atom.kind = self.connectivity
        out_atom.kind = self.kind
        return out_atom

    def translate(self, x1, y1, z1):
        """Translate the atom by some vector."""
        self.x += x1
        self.y += y1
        self.z += z1
        return

    def v_translate(self, vec_trans):
        """Translate the atom by some vector."""
        self.x += vec_trans[0]
        self.y += vec_trans[1]
        self.z += vec_trans[2]
        return

    def set_connectivity(self, in_atoms, in_row):
        """
        Set the connectivity and the kind of the atom.

        This function needs a row of a connectivity matrix which can be obtained
        with functions from assign_charges.py

        Check the constructor at the top of this file for more info on connectivity
        and kind.

        Parameters
        ----------
        in_atoms : list of atoms
            Atoms in the system of which this atom is a part
        in_row : 1-d array-like
            The row of the connectivity matrix of in_atoms which corresponds to
            this atom

        """
        links = []
        for i, atom in enumerate(in_atoms):
            if in_row[i] != 0:
                links.append((in_atoms[i].elem, in_row[i]))
        self.connectivity = frozenset(Counter(links).most_common())
        self.kind = (self.elem, self.connectivity)
        return

    def es_pot(self, position):
        """
        Return the electorstatic potential generated by this Atom

        Parameters
        ----------
        position : 3x1 np array
            The point at which the potential should be evaluated
        Returns
        -------
        pot : float
            The potential of the Atom felt at the input position

        """
        pot = self.q / self.v_dist(position)
        return pot
