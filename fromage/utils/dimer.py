"""
The class represents pairs of Mol objects
"""
import fromage.utils.array_operations as ao
import numpy as np
from scipy.spatial.distance import cdist

def make_dimer(mol_a,mol_b):
    """
    Build a Dimer object

    Parameters
    ----------
    mol_a,mol_b : Mol objects
        The molecules constituting the dimer
    Returns
    -------
    dim_out : Dimer object
        The dimer of the two molecules

    """
    dim_out = Dimer(mol_a,mol_b)
    return dim_out

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

    def write_xyz(self, name):
        """Write an xyz file of the Dimer"""
        whole_mol = self.mol_a + self.mol_b
        whole_mol.write_xyz(name)

    def identical_to(self,other_dimer):
        """Check if the dimer is the same as another one"""
        result = False
        if self.mol_a.same_atoms_as(other_dimer.mol_a):
            if self.mol_b.same_atoms_as(other_dimer.mol_b):
                result = True
        elif self.mol_a.same_atoms_as(other_dimer.mol_b):
            if self.mol_b.same_atoms_as(other_dimer.mol_a):
                result = True
        return result

    def angles(self):
        """
        Return the three descriptor angles of the dimer

        Returns
        -------
        out_arr : 3 x 1 numpy array
            The three angles alpha, beta, gamma

        """
        if np.count_nonzero(self.mol_a.geom.perp_ax) == 0:
            self.mol_a.calc_axes()
        if np.count_nonzero(self.mol_b.geom.perp_ax) == 0:
            self.mol_b.calc_axes()
        out_lis = [ao.vec_angle(self.mol_a.geom.prin_ax,self.mol_b.geom.prin_ax),
                    ao.vec_angle(self.mol_a.geom.sec_ax,self.mol_b.geom.sec_ax),
                    ao.vec_angle(self.mol_a.geom.perp_ax,self.mol_b.geom.perp_ax)]
        out_arr = np.array(out_lis)

        return out_arr

    def slip_angle(self):
        """
        Return the slip angles for the dimer

        If the dimer is arranged face-to-face, their centroids can be said
        to be more or less slipped (deviating from being aligned). A measure
        of this slip is the angle between the normal axis of one of the monomers
        and the centroid-centroid vector. Note that the slip angle of dimer IJ
        is not related to that of JI so both are returned. Depending on the
        direction of the vectors, two definitions can be chosen for the angle
        (both of them adding up to 180). We return the smalles one of the two.

        Returns
        -------
        slip_angle : floats
            The slip angle

        """
        if np.count_nonzero(self.mol_a.geom.perp_ax) == 0:
            self.mol_a.calc_axes()
        if np.count_nonzero(self.mol_b.geom.perp_ax) == 0:
            self.mol_b.calc_axes()
        # vector gonig from A to B
        cen_cen = self.mol_b.centroid() - self.mol_a.centroid()
        # test angle with positive cen_cen and negative cen_cen. keep smallest
        slip_angle_a_pos = ao.vec_angle(self.mol_a.geom.perp_ax, cen_cen)
        slip_angle_a_neg = 180 - slip_angle_a_pos
        slip_angle_b_pos = ao.vec_angle(self.mol_b.geom.perp_ax, cen_cen)
        slip_angle_b_neg = 180 - slip_angle_b_pos
        slip_angle_a = min((slip_angle_a_pos,slip_angle_a_neg))
        slip_angle_b = min((slip_angle_b_pos,slip_angle_b_neg))

        final_slip = min(slip_angle_a, slip_angle_b)
        return final_slip

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
        Return the 27 images of the dimer produced by translating monomer 2

        Monomer 2 is translated in every combination of -a,0,a ; -b,0,b and
        -c,0,c. This includes the (0,0,0) translation which is the original
        dimer.

        Parameters
        ----------
        vectors : 3 x 3 numpy array
            Vectors of the lattice periodicity
        Returns
        -------
        images : list of 27 Dimer objects
            The 27 images of the dimer, always with mol_a remaining in first
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

    def sorted_inter_distances(self):
        """
        Return the sorted intermolecular atomic distances

        That is, only distances between an atom belonging to mol_a and an atom
        belonging to mol_b

        Returns
        -------
        inter_distances : Nat x Nat numpy array
            The sorted distance matrix of interatomic distances across molecules

        """
        arr_a = self.mol_a.coord_array()
        arr_b = self.mol_b.coord_array()

        unsorted = cdist(arr_a,arr_b)
        unsort_arr = unsorted.flatten()
        inter_distances = np.sort(unsort_arr)

        return inter_distances

    def same_geom(self,other,tol=10e-4):
        """
        Checks whether the dimer has the same geometry as another

        Atom types are ignored. The comparison is done by sorting all inter
        molecular distances (atom to atom). This fingerprint should be
        sufficient in general. The rmsd between dimers is checked and needs
        to fall below a threshold

        Parameters
        ----------
        other : Mol object
            Molecule to compare with
        tol : float
            Tolerance for the RMSD of the two sorted distance arrays
        Returns
        -------
        same : bool
            True if the molecules have the same geometry, else False

        """
        dists_a = self.sorted_inter_distances()
        dists_b = other.sorted_inter_distances()

        same = ao.rmsd(dists_a,dists_b) < tol

        return same

    def mols_are_linear(self):
        """State that the geometries of the molecules are linear-like"""
        self.mol_a.geom.linear = True
        self.mol_b.geom.linear = True
