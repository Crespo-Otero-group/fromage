"""Functions relating to numpy arrays"""
import numpy as np
from scipy.spatial.distance import cdist

from ._planes import plane_from_coord, quadrangle_from_coord, embedded_vert, project_point, project_pair_to_vector, project_quad_to_vectors
from ._matrix import cross_product_matrix, rotation_matrix

def distance(vector_1, vector_2):
    """Return the distance between two points"""
    dis = np.linalg.norm(vector_1-vector_2)
    return dis

def vec_angle(vector_1, vector_2, degrees = True):
    """
    Return the angle between two numpy vectors.

    The angle is in the range [0,180] for degrees or [0,1] for
    radians. arctan is used instead of e.g. arccos as it has a more robust
    definition for edge cases where floating point numbers might get
    arccos(1.000001) = NaN

    Parameters
    ----------
    vector_1, vector_2 : N x 1 numpy array
        The vectors whose angle needs to be calculated
    degrees : bool (optional)
        Result in degrees or in radians. Default = True, so degrees
    Returns
    -------
    out_angle : float
        The angle between vectors
    """
    dot = np.dot(vector_1,vector_2)
    cross_norm = np.linalg.norm(np.cross(vector_1,vector_2))
    ang = np.arctan2(cross_norm,dot)
    if degrees:
        ang = np.degrees(ang)
    return ang

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

def rmsd(array_a, array_b):
    """
    Calculate the RMSD between two 1d arrays

    Parameters
    ----------
    array_a, array_b : 1d numpy arrays
        The arrays to be compared
    Returns
    -------
    rmsd : float
        The Root Mean Square Deviation of the elements of the array

    """
    diff = array_a - array_b
    diff2 = np.square(diff)
    diff2_sum = np.sum(diff2)
    norm_diff2_sum = diff2_sum/len(array_a)
    rmsd = np.sqrt(norm_diff2_sum)

    return rmsd

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

def orthogonalise_sym(vectors):
    """
    Return two orthogonal vectors based on the original ones, both rotated

    We wish to orthogonalise two vectors but moving each one by the same amount.
    That is the angle of rotation of vec_1 should be negative that of vec_2 and
    the two ending vectors should be orthogonal.

    Parameters
    ----------
    vectors : 2 x 3 numpy array
        The vectors to be orthogonalised
    Returns
    -------
    o_vecs: 2 x 3 numpy array
        Orthogonalised vectors
    """
    ang = vec_angle(vectors[0],vectors[1])
    remainder = 90 - ang
    disp = remainder/2
    perp_unnormal = np.cross(vectors[0],vectors[1])
    normal = perp_unnormal / np.linalg.norm(perp_unnormal)

    rot_1 = rotation_matrix(normal,-disp)
    rot_2 = rotation_matrix(normal,disp)

    ovec_1 = np.dot(rot_1,vectors[0])
    ovec_2 = np.dot(rot_2,vectors[1])

    o_vecs = np.array([ovec_1,ovec_2])
    return o_vecs

def orthogonalise_asym(vectors):
    """
    Return two orthogonal vectors based on the original ones, both rotated

    Here, we only rotate the second vector in order to make it perpendicular to
    the first one

    Parameters
    ----------
    vectors : 2 x 3 numpy array
        The vectors to be orthogonalised
    Returns
    -------
    o_vecs: 2 x 3 numpy array
        Orthogonalised vectors where the second one has been rotated to become
        orthogonal
    """
    ang = vec_angle(vectors[0],vectors[1])
    remainder = 90 - ang
    disp = remainder
    perp_unnormal = np.cross(vectors[0],vectors[1])
    normal = perp_unnormal / np.linalg.norm(perp_unnormal)

    rot_2 = rotation_matrix(normal,disp)

    ovec_1 = vectors[0]
    ovec_2 = np.dot(rot_2,vectors[1])

    o_vecs = np.array([ovec_1,ovec_2])
    return o_vecs
