import numpy as np
import fromage.utils.array_operations as ao

def coord_array(self):
    """
    Return a numpy array of the coordinates

    Returns
    -------
    coord_arr : Nat x 3 numpy array
        Array of the form [[x1,y1,z1],[x2,y2,z2],...]

    """
    nat = len(self)
    coord_arr = np.zeros((nat,3))
    for i,atom in enumerate(self):
        coord_arr[i][0] = atom.x
        coord_arr[i][1] = atom.y
        coord_arr[i][2] = atom.z
    return coord_arr

def extreme_at_pairs(self, n_pairs):
    """
    Return a list of pairs of extreme atom coordinates

    Parameters
    ----------
    n_pairs : int
        Number of extreme atom pairs requested
    Returns
    -------
    pairs : numpy array of n_pairs x 2 x 3

    """
    coords = self.coord_array()
    dist_mat = ao.dist_mat(coords)
    pairs_inds = ao.find_largest(dist_mat,2)
    pairs = np.zeros((n_pairs,2,3))
    for i,ind in enumerate(pairs_inds):
        pairs[i][0] = coords[ind[0]]
        pairs[i][1] = coords[ind[1]]

    return pairs
