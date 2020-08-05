"""Functions relating to numpy arrays"""
import numpy as np
import itertools
from scipy.spatial.distance import cdist
from fromage.utils.atom import Atom

from ._planes import (
    plane_from_coord,
    quadrangle_from_coord,
    embedded_vert,
    project_point,
    project_pair_to_vector,
    project_quad_to_vectors,
)
from ._matrix import cross_product_matrix, rotation_matrix, reflection_matrix


def distance(vector_1, vector_2):
    """Return the distance between two points"""
    dis = np.linalg.norm(vector_1 - vector_2)
    return dis


def vec_angle(vector_1, vector_2, degrees=True):
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
    dot = np.dot(vector_1, vector_2)
    cross_norm = np.linalg.norm(np.cross(vector_1, vector_2))
    ang = np.arctan2(cross_norm, dot)
    if degrees:
        ang = np.degrees(ang)
    return ang


def closest(reference, points):
    """Return the closest point to the reference and the distance to it"""
    min_dis = float("inf")
    for point in points:
        dis = distance(reference, point)
        if dis < min_dis:
            min_dis = dis
            closest_point = point
    return closest_point, min_dis


def furthest(reference, points):
    """Return the furthest point to the reference and the distance to it"""
    max_dis = -float("inf")
    for point in points:
        dis = distance(reference, point)
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
    dist_mat = np.tril(cdist(in_array, in_array))
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
    rmsd_val : float
        The Root Mean Square Deviation of the elements of the array

    """
    diff = array_a - array_b
    diff2 = np.square(diff)
    diff2_sum = np.sum(diff2)
    norm_diff2_sum = diff2_sum / len(array_a)
    rmsd_val = np.sqrt(norm_diff2_sum)

    return rmsd_val


def coord_rmsd(array_a, array_b):
    """
    Calculate the distance RMSD between sets of coordinates

    Parameters
    ----------
    array_a, array_b : N x 3 numpy arrays
        The arrays [[x1, y1, z1], [x2, y2, z2], ...]
    Returns
    -------
    rmsd_val : float
        The Root Mean Square Deviation of the distance between points

    """
    displacement_vectors = array_a - array_b
    distances = np.linalg.norm(displacement_vectors, axis=1)
    # rmsd with 0 distances
    rmsd_val = rmsd(np.zeros(len(distances)), distances)

    return rmsd_val

    return rmsd_val


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
        folded_index = np.unravel_index(flat_index, shape)
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
    ang = vec_angle(vectors[0], vectors[1])
    remainder = 90 - ang
    disp = remainder / 2
    perp_unnormal = np.cross(vectors[0], vectors[1])
    normal = perp_unnormal / np.linalg.norm(perp_unnormal)

    rot_1 = rotation_matrix(normal, -disp)
    rot_2 = rotation_matrix(normal, disp)

    ovec_1 = np.dot(rot_1, vectors[0])
    ovec_2 = np.dot(rot_2, vectors[1])

    o_vecs = np.array([ovec_1, ovec_2])
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
    ang = vec_angle(vectors[0], vectors[1])
    remainder = 90 - ang
    disp = remainder
    perp_unnormal = np.cross(vectors[0], vectors[1])
    normal = perp_unnormal / np.linalg.norm(perp_unnormal)

    rot_2 = rotation_matrix(normal, disp)

    ovec_1 = vectors[0]
    ovec_2 = np.dot(rot_2, vectors[1])

    o_vecs = np.array([ovec_1, ovec_2])
    return o_vecs


def array2atom(template, pos):
    """
    Turn an array of the form x1, y1, z1, x2, y2, z2 etc. into a list of Atom
    objects

    Parameters
    ----------
    template : list of Atom objects
        A list of the same length of the desired one used to determine the
        elements of the atoms
    pos : list of floats
        List of coordinates of the form x1, y1, z1, x2, y2, z2 etc.
    Returns
    -------
    out_atoms : list of Atom objects
        Resulting atoms

    """
    sliced_pos = [pos[i : i + 3] for i in range(0, len(pos), 3)]
    out_atoms = []
    for atom in zip(template, sliced_pos):
        new_atom = Atom(atom[0].elem, atom[1][0], atom[1][1], atom[1][2], 0)
        out_atoms.append(new_atom)
    return out_atoms


def possible_translations(lat_vectors):
    """
    Return possible first order translations from a set of lattice vectors

    Includes possitive and negative translations

    Parameters
    ----------
    lat_vectors : 3 x 3 numpy array
        Lattice vectors
    Returns
    -------
    possible_trans : 27 x 3 numpy array
        The 27 translations that can be made by adding together lattice vectors

    """
    multi = np.array([-1, 0, 1])
    # NB the iterable needs to become a list before an array because the iterable
    # loses its elements as soon as they get accessed, and arrays parse the
    # elements by accessing them several times
    multi_sets = np.array(list(itertools.product(multi, multi, multi)))
    possible_trans = np.einsum("ij,jk->ij", multi_sets, lat_vectors)

    return possible_trans


def dist_vec(coord_a, coord_b):
    """
    Return the vectors associate with the distance matrix between two sets of coordinates

    The vectors are in the direction A->B

    Parameters
    ----------
    coord_a, coord_b : (M or N) x 3 numpy arrays
        The coordinate arrays [[x1, y1, z1], [x2, y2, z2], ...]

    Returns
    -------
    displacements : M x N x 3 numpy array
        The matrix containing the vectors between all points

    """
    displacements = coord_b[np.newaxis, :, :] - coord_a[:, np.newaxis, :]

    return displacements


def per_dist_mat(coord_a, coord_b, lat_vec=None):
    """
    Return the distance matrix between two sets of coordinates

    There is an option to do so within a periodic cell.

    Parameters
    ----------
    coord_a, coord_b : (M or N) x 3 numpy arrays
        The coordinate arrays [[x1, y1, z1], [x2, y2, z2], ...]
    lat_vec: 3 x 3 numpy array or None
        If lattice vectors are supplied, the distance vectors are chosen such that the
        distance is the lowest, considering periodic images

    Returns
    -------
    dis_mat : M x N numpy array
        The distance matrix

    """
    displacements = dist_vec(coord_a, coord_b)

    # not periodic case
    if lat_vec is None:
        dis_mat = np.linalg.norm(displacements, axis=-1)
    # periodic case
    else:
        # get an array where, for each displacement, an additional axis of size
        # 27 gives its results with the added possible translations
        possible_disp = displacements[:, :, np.newaxis, :] + possible_translations(
            lat_vec
        )
        # get the norms, in the shape M x N x 27
        possible_norms = np.linalg.norm(possible_disp, axis=-1)
        # pick the smallest of the 27 for each distance
        dis_mat = np.min(possible_norms, axis=-1)

    return dis_mat


def dir_to_frac(coords, lat_vec):
    """
    Return coordinates in fractional coordinates

    Parameters
    ----------
     coords : N x 3 numpy array
        Direct coordinate array
    lat_vec: 3 x 3 numpy array or None
        Lattice vectors
    Returns
    -------
    frac_coords : N x 3 numpy array
        Fractional coordinate array

    """
    # transpose lattice vectors
    M = np.transpose(lat_vec)
    # inverse transformation matrix
    U = np.linalg.inv(M)
    # apply transform and modulo 1
    frac_coords = np.mod(np.dot(U, coords.T).T, 1)

    return frac_coords


def frac_to_dir(coords, lat_vec):
    """
    Return coordinates in direct coordinates

    Parameters
    ----------
     coords : N x 3 numpy array
        Fractional coordinate array
    lat_vec: 3 x 3 numpy array or None
        Lattice vectors
    Returns
    -------
    dir_coords : N x 3 numpy array
        Direct coordinate array

    """
    dir_coords = np.matmul(lat_vec.T, coords.T).T

    return dir_coords


def confine(coords, lat_vec):
    """
    Translate all coordinates so that they fit in the conventional unit cell

    Parameters
    ----------
    coords : N x 3 numpy array
        Input coordinates
    lat_vec : 3 x 3 numpy array
        Lattice vectors
    Returns
    -------
    confined_coords = N x 3 numpy array
        Coordinates in direct space after being translated in the cell

    """
    confined_coords = frac_to_dir(dir_to_frac(coords, lat_vec), lat_vec)

    return confined_coords


def per_translate(coords, trans, lat_vec):
    """
    Translate atoms and then confine them to a conventional cell

    Parameters
    ----------
    coords : N x 3 numpy array
        Coordinate array
    trans : length 3 numpy array
        Translation vector
    lat_vec: 3 x 3 numpy array or None
        Lattice vectors
    Returns
    -------
    new_coords : N x 3 numpy array
        Coordinate array after translation

    """
    new_coords = confine(coords + trans, lat_vec)

    return new_coords
