"""
The class represents pairs of Mol objects
"""
import fromage.utils.array_operations as ao
import numpy as np
class Dimer(object):
    """
    Object representing a pair of molecules

    Attributes
    ----------
    mols : list of two Mol objects
        The two molecules constituting the dimer
    alpha, beta, gamma : floats
        The describing angles. alpha is the angle between principal axes,
        beta secondary and gamma perpendicular

    """

    def __init__(self, mol_a=None, mol_b=None):
        self.mol_a = mol_a
        self.mol_b = mol_b
        self.alpha = None
        self.beta = None
        self.gamma = None

    def __repr__(self):
        out_str = "Mol A\n" + self.mol_a.__str__() + "Mol B\n" + self.mol_b.__str__()
        return out_str

    def __str__(self):
        return self.__repr__()

    def angles(self):
        """
        Return the three descriptor angles of the dimer

        Returns
        -------
        out_arr : 3 x 1 numpy array
            The three angles alpha, beta, gamma

        """
        if self.mol_a.geom.perp_ax == None:
            self.mol_a.calc_axes()
        if self.mol_b.geom.perp_ax == None:
            self.mol_b.calc_axes()
        out_lis = [ao.vec_angle(self.mol_a.geom.prin_ax,self.mol_b.geom.prin_ax),
                    ao.vec_angle(self.mol_a.geom.sec_ax,self.mol_b.geom.sec_ax),
                    ao.vec_angle(self.mol_a.geom.perp_ax,self.mol_b.geom.perp_ax)]
        out_arr = np.array(out_lis)

        return out_arr

    def calc_angles(self):
        """Set the three descriptor angles"""
        descriptor_angles = self.angles()
        self.alpha = descriptor_angles[0]
        self.beta = descriptor_angles[1]
        self.gamma = descriptor_angles[2]

        return

    def inter_distance(self, method='atomic',mode='dis'):
        """
        Return the intermolecular distance

        This can be defined as centroid-centroid or closest atom-atom. The
        latter can further be refined by taking into account van der waals
        or covalent radii.

        Parameters
        ----------
        method : str (optional)
            'atomic' or 'centroid' respectively correspond to closest atom-atom
            distance and centroid-centroid distance. Default 'atomic'
        mode : str (optional)
            'dis', 'cov or 'vdw' determines the scaling when the method is
            atomic.
        Returns
        -------
        dis : float
            The intermolecular distance

        """
        if method == 'centroid':
            dis = self.inter_dist_centroid()
        elif method == 'atomic':
            dis = self.inter_dist_atomic(mode=mode)
        else:
            raise ValueError("The only methods available are 'atomic' and 'centroid'.\
        You requested: " + str(method))

        return dis

    def inter_dist_centroid(self):
        """Return distance between centroids of constituent fragments"""
        cen_1 = self.mol_a.centroid()
        cen_2 = self.mol_b.centroid()

        diff = cen_1 - cen_2

        dis = np.linalg.norm(diff)

        return dis

    def inter_dist_atomic(self, mode='dis'):
        """
        Return the closest atomic distance between molecules

        Parameters
        ----------
        mode : str (optional)
            'dis', 'cov or 'vdw' determines the scaling when the method is
            atomic.
        Returns
        -------
        dis : float
            The intermolecular distance, possibly tweaked by covalent or vdw
            radii

        """
        min_dist = float('inf')

        for atom_i in self.mol_a:
            for atom_j in self.mol_b:
                tmp_dis = atom_i.dist(atom_j,ref=mode)
                if tmp_dis < min_dist:
                    min_dist = tmp_dis

        return min_dist

    def images(self, vectors):
        """
        Return the 9 images of the dimer produced by translating monomer 2

        Monomer 2 is translated in every combination of -a,0,a ; -b,0,b and
        -c,0,c. This includes the (0,0,0) translation which is the original
        dimer.

        Parameters
        ----------
        vectors : 3 x 3 numpy array
            Vectors of the lattice periodicity
        Returns
        -------
        images : list of 9 Dimer objects
            The 9 images of the dimer, always with mol_a remaining in first
            position.

        """
        translations = [-1,0,1]
        static_mol = self.mol_a
        moving_mol = self.mol_b
        images = []
        for tra_a in translations:
            for tra_b in translations:
                    for tra_c in translations:
                        vector = vectors[0] * tra_a + \
                                 vectors[1] * tra_b + \
                                 vectors[2] * tra_c
                        mol_image = moving_mol.translated(vector)
                        new_dimer = Dimer(static_mol,mol_image)
                        images.append(new_dimer)
        return images
