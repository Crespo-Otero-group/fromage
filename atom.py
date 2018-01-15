"""Defines the Atom object"""

import numpy as np
import periodic as per
from collections import Counter


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

    """

    def __init__(self, elemIn="H", xIn=0.0, yIn=0.0, zIn=0.0, qIn=0.0, num=1):
        self.elem = elemIn
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.q = 0.0
        self.num = 1
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
        self.vdw = table[self.elem.lower()]["vdw"]


        # to string methods to be used mainly for debugging and .qc file
    def __repr__(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(self.elem, self.x, self.y, self.z, self.q)

    def __str__(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(self.elem, self.x, self.y, self.z, self.q)

        # equality function
    def __eq__(self, other):
        return self.elem.lower() == other.elem.lower() and self.x == other.x and self.y == other.y and self.z == other.z and self.q == other.q

    def very_close(self, other):
        """Check if two atoms are very close together"""
        thresh = 0.001
        x_cond = abs(self.x - other.x) < thresh
        y_cond = abs(self.y - other.y) < thresh
        z_cond = abs(self.z - other.z) < thresh

        cond = x_cond and y_cond and z_cond
        return cond

    def xyz_str(self):
        """Return a string of the atom in xyz format"""
        return "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(self.elem, self.x, self.y, self.z)

    def dist2(self, x1, y1, z1):
        """Return distance squared of the atom from a point"""
        r = (self.x - x1) ** 2 + (self.y - y1) ** 2 + (self.z - z1) ** 2
        return r

    def dist(self, x1, y1, z1):
        """Return distance of the atom from a point"""
        r = np.sqrt(self.dist2(x1, y1, z1))
        return r

    def dist_lat(self, x1, y1, z1, aVec, bVec, cVec):
        """
        Find the shortest distance to a point in a periodic system.

        Parameters
        ----------
        x1,y1,z1 : floats
            Cartesian coordinates of the target point
        aVec,bVec,cVec : 3x1 array-likes
            Unit cell vectors

        Returns
        -------
        rMin : float
            Minimal distance to the point
        x3,y3,z3 : floats
            Coordinates of the closest image to the point

        """
        # null vector
        nVec = (0, 0, 0)
        # negative vectors
        aVecN = [-i for i in aVec]
        bVecN = [-i for i in bVec]
        cVecN = [-i for i in cVec]

        # sets comprised of the lattice vector,
        # the null vector and the negative lattice vector
        aSet = [aVec, nVec, aVecN]
        bSet = [bVec, nVec, bVecN]
        cSet = [cVec, nVec, cVecN]

        # minimum r distance
        rMin = float("inf")

        # loop over all possible translations of the input point
        for trans1 in aSet:
            for trans2 in bSet:
                for trans3 in cSet:
                    x2 = x1 + trans1[0] + trans2[0] + trans3[0]
                    y2 = y1 + trans1[1] + trans2[1] + trans3[1]
                    z2 = z1 + trans1[2] + trans2[2] + trans3[2]
                    r = np.sqrt((self.x - x2) ** 2 + (self.y - y2)
                                ** 2 + (self.z - z2) ** 2)
                    # if this particular translation of the point is the closest
                    # to the atom so far
                    if r < rMin:
                        rMin = r
                        # image coordinates
                        x3 = x2
                        y3 = y2
                        z3 = z2
        return rMin, x3, y3, z3

    def translated(self, x1, y1, z1):
        "Return a new atom which is a translated copy."
        xout, yout, zout = self.x, self.y, self.z
        xout += x1
        yout += y1
        zout += z1
        outAtom = Atom(self.elem, xout, yout, zout, self.q)
        return outAtom

    def translate(self, x1, y1, z1):
        "Translate the atom by some vector."
        self.x += x1
        self.y += y1
        self.z += z1
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
