"""Functions relating to numpy arrays"""
import numpy as np
from scipy.spatial.distance import cdist

from ._planes import plane_from_coord, quadrangle_from_coord, embedded_vert, project_point, project_pair_to_vector, project_quad_to_vectors
from ._angles import vec_angle

def distance(vector_1, vector_2):
    """Return the distance between two points"""
    dis = np.linalg.norm(vector_1-vector_2)
    return dis

def closest(reference,points):
    """Return the closest point to the reference and the distance to it"""
    min_dis = float('inf')
    for point in points:
        dis = distance(reference,point)
        if dis < min_dis:
            min_dis = dis
            closest_point = point
    return closest_point, min_dis

def furthest(reference,points):
    """Return the furthest point to the reference and the distance to it"""
    max_dis = -float('inf')
    for point in points:
        dis = distance(reference,point)
        if dis > max_dis:
            max_dis = dis
            closest_point = point
    return closest_point, max_dis

def dist_mat(in_array):
    """
    Return lower triangular distance matrix of coordinate array

    Parameters
    ----------
    in_array : Nat X 3 numpy array
        Coordinate array
    Returns
    -------
    dist_mat : Nat X Nat numpy array
        Lower triangular distance matrix

    """
    dist_mat = np.tril(cdist(in_array,in_array))
    return dist_mat

def find_largest(in_array, n_largest):
    """
    Return the coordinates of the N largest elements of an ND-array

    Parameters
    ----------
    in_array : M x N numpy array
        Numpy array
    n_largest : int
        Number of required pairs of atoms
    Returns
    -------
    n_largest : n_largest X 2 array-like
        The coordinates of the n_largest elements. It's an n_largest dimensional
        array of tuples.

    """
    new_arr = np.copy(in_array)
    shape = np.shape(new_arr)
    indices = []
    while len(indices) < n_largest:
        flat_index = np.argmax(new_arr)
        folded_index = np.unravel_index(flat_index,shape)
        indices.append(folded_index)
        new_arr[folded_index] = 0
    return indices
