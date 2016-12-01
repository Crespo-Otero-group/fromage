# Atom class to interface between formatting
# in different programs
from numpy import *


class Atom:
        # Atom class has element type, coordinates and charge

    def __init__(self, elemIn="H", xIn=0.0, yIn=0.0, zIn=0.0, qIn=0.0):
        self.elem = elemIn
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.q = 0.0
        try:
            self.x = float(xIn)
            self.y = float(yIn)
            self.z = float(zIn)
            self.q = float(qIn)

        except ValueError:
            print "Some coordinates or charges cannot be cast to float!"

        # to string methods to be used mainly for debugging and .qc file
    def __repr__(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(self. elem, self.x, self.y, self.z, self.q)
        # return str(self.elem) + "\t" + str(self.x) + "\t" + str(self.y) +
        # "\t" + str(self.z) + "\t" + str(self.q) + "\n"

    def __str__(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(self. elem, self.x, self.y, self.z, self.q)
        # return str(self.elem) + "\t" + str(self.x) + "\t" + str(self.y) +
        # "\t" + str(self.z) + "\t" + str(self.q) + "\n"

        # writes the atom coordinates in xyz format

    def xyzStr(self):
        return "{:>6} {:10.6f} {:10.6f} {:10.6f}".format(self. elem, self.x, self.y, self.z)
        # return str(self.elem) + "\t" + str(self.x) + "\t" + str(self.y) +
        # "\t" + str(self.z) + "\n"

        # equality function
    def __eq__(self, other):
        return self.elem == other.elem and self.x == other.x and self.y == other.y and self.z == other.z and self.q == other.q

        # returns the distance of the atom from an input point
    def dist(self, x1, y1, z1):
        r = sqrt((self.x - x1) ** 2 + (self.y - y1) ** 2 + (self.z - z1) ** 2)
        return r

        # returns the distance of the atom from an input point or
        # its closest periodic image and the coordinates of
        # whichever was closest
    def distLat(self, x1, y1, z1, aVec, bVec, cVec):
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
                    r = sqrt((self.x - x2) ** 2 + (self.y - y2)
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

        # returns an atom translated by some vector
    def translate(self, x1, y1, z1):
        xout, yout, zout = self.x, self.y, self.z
        xout += x1
        yout += y1
        zout += z1
        outAtom = Atom(self.elem, xout, yout, zout, self.q)
        return outAtom

    def electrons(self):
        total=0
        valence=0

        element = self.elem.lower()
        if element == "h":
            total = 1
            valence = 1
        elif element == "c":
            total = 6
            valence = 4
        elif element == "n":
            total = 7
            valence = 5
        elif element == "o":
            total = 8
            valence = 6

        return (valence,total)
