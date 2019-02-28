"""Functions pertaining to numpy arrays. Many functions operate on coordinate
arrays of the form [[x1,y1,z1],[x2,z2,y2],...] and are later on used in
mol._geom"""
import numpy as np
from scipy.spatial.distance import cdist

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

def plane_from_coord(coord_arr):
    """
    Return plane coefficients which best includes the coordinates

    This is done via a singular value decomposition of the coordinates. Read
    the following:
    http://caves.org/section/commelect/DUSI/openmag/pdf/SphereFitting.pdf

    Parameters
    ----------
    coord_arr : Nat x 3 numpy array
        The input coordinates array
    Returns
    -------
    a, b, c, d : floats
        Plane equation coefficients such that a point on the plane is:
        ax + by + cz + d = 0

    """
    centroid = np.average(coord_arr, axis=0)
    # coordinates translated to origin
    cen_coords = coord_arr - centroid

    U, S, Vt = np.linalg.svd(cen_coords)
    # The last row of V matrix indicate the eigenvectors of
    # smallest eigenvalues (singular values).
    N = Vt[-1]

    # Extract a, b, c, d coefficients.
    x0, y0, z0 = centroid
    a, b, c = N
    d = -(a * x0 + b * y0 + c * z0)

    return a, b, c, d

def extreme_pairs(coord_arr, n_pairs):
    """
    Return a list of pairs of extreme coordinates

    Parameters
    ----------
    coord_arr : Nat x 3 numpy array
        The input coordinates array
    n_pairs : int
        Number of extreme atom pairs requested
    Returns
    -------
    pairs : numpy array of n_pairs x 2 x 3
        Coordinates of the N pairs of points

    """
    dmat = dist_mat(coord_arr)
    pairs_inds = find_largest(dmat,2)
    pairs = np.zeros((n_pairs,2,3))
    for i,ind in enumerate(pairs_inds):
        pairs[i][0] = coord_arr[ind[0]]
        pairs[i][1] = coord_arr[ind[1]]

    return pairs

def embedded_pairs(coord_pairs):
    """
    Return the atom pairs defining the quadrilateral embedded in the input

    The input is two pairs of atoms where each pair is the diagonal of the
    quadrilateral, defined as the two most distant coordinate pairs. The output
    is the corresponding pair.

    Parameters
    ----------
    coord_pairs : numpy array of 2 x 2 x 3
        Pairs of coordinates
    Returns
    -------
    out_pairs : numpy array of 2 x 2 x 3
        Pairs of coordinates where the first pair is the longer diagonal of the
        new quadrilateral

    """
    first_pair = coord_pairs[0]
    second_pair = coord_pairs[1]

    # define the coordinates of the diagonal A
    extremum_a_1 = np.mean([first_pair[0],second_pair[0]],axis=0)
    extremum_a_2 = np.mean([first_pair[1],second_pair[1]],axis=0)
    # same for the diagonal B
    extremum_b_1 = np.mean([first_pair[0],second_pair[1]],axis=0)
    extremum_b_2 = np.mean([first_pair[1],second_pair[0]],axis=0)

    # now determine which diagonal is the long one

    # measure the distance of atom 1 of pair 1 to atom 1 of pair 2
    dis_1 = np.linalg.norm(first_pair[0] - second_pair[0])
    # again for atom 1 of pair to atom 2 of pair 2
    dis_2 = np.linalg.norm(first_pair[0] - second_pair[1])

    # if firstpair_[0], second_pair[0] is the short side of the original
    # quadrilateral
    if dis_1 < dis_2:
        long_1 = extremum_a_1
        long_2 = extremum_a_2
        short_1 = extremum_b_1
        short_2 = extremum_b_2
    # if it's the other way rond
    else:
        long_1 = extremum_b_1
        long_2 = extremum_b_2
        short_1 = extremum_a_1
        short_2 = extremum_a_2

    out_pairs = np.array([[long_1, long_2],[short_1, short_2]])

    return out_pairs
