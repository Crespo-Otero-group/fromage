"""The class represents pairs of Mol objects
"""

class Dimer(object):
    """
    Object representing a pair of molecules

    Attributes
    ----------
    mols : list of two Mol objects
        The two molecules constituting the dimer
    """

    def __init__(self, mols=[]):
        self.mols = mols
    def __repr__(self):
        out_str = "Mol A\n" + self.mols[0].__str__() + "Mol B\n" + self.mols[1].__str__()
        return out_str
    def __str__(self):
        return self.__repr__()
