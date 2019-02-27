import numpy as np
from scipy.spatial.distance import cdist

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

def pairwise_distances(self):
    """
    Return a numpy array of the pairwise squared distances

    Returns
    -------
    dis_arr : Nat x Nat numpy array
        Squared distance matrix (lower triangular and 0 diagonal)

    """
    coords = self.coord_array()
    dis_arr = cdist(coords,coords)

    return dis_arr
