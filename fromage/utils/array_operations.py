"""Functions pertaining to numpy arrays. Many functions operate on coordinate
arrays of the form [[x1,y1,z1],[x2,z2,y2],...] and are later on used in
mol._geom"""
import numpy as np
from scipy.spatial.distance import cdist

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

def plane_from_coord(coord_arr):
    """
    Return plane coefficients which best includes the coordinates

    This is done via a singular value decomposition of the coordinates. Read
    the following:
    http://caves.org/section/commelect/DUSI/openmag/pdf/SphereFitting.pdf
    Given the plane equation ax + by + cz + d = 0, this can be interpeted as
    a unit vector (a,b,c), perpendicular to the plane and an origin-plane
    distance of d

    Parameters
    ----------
    coord_arr : Nat x 3 numpy array
        The input coordinates array
    Returns
    -------
    res : length 3 numpy array
        Plane equation coefficients such that a point on the plane is:
        ax + by + cz + d = 0. The array is [a,b,c,d]

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

    res = np.array([a,b,c,d])
    return res

def quadrangle_from_coord(coord_arr):
    """
    Return a list of ordered vertices of a quadrangle from extreme points

    The quadrangle has, for diagonals, the two largest inter - coordinate distances.
    Given the quadrangle ABCD the output is [A, B, C, D] such that
    AC > BD > any other distance in the coordinate array. Furthermore, AB is the
    longest side. This completely determines the quadrilateral.

    Parameters
    ----------
    coord_arr : Nat x 3 numpy array
        The input coordinates array
    Returns
    -------
    vertices_out : 4 x 1 numpy array
        Coordinates of the four extreme atoms

    """
    dmat = dist_mat(coord_arr)
    pairs_inds = find_largest(dmat,2)
    # the first two are either A or C and the last two B or D
    # We label [AC1,AC2,BD1,BD2]
    unordered = [coord_arr[pairs_inds[0][0]],
                coord_arr[pairs_inds[0][1]],
                coord_arr[pairs_inds[1][0]],
                coord_arr[pairs_inds[1][1]]]

    # calculate the four sides
    AC1_BD1_dis = distance(unordered[0],unordered[2])
    AC1_BD2_dis = distance(unordered[0],unordered[3])
    AC2_BD1_dis = distance(unordered[1],unordered[2])
    AC2_BD2_dis = distance(unordered[1],unordered[3])

    # find the longest side
    max_side = max([AC1_BD1_dis,AC1_BD2_dis,AC2_BD1_dis,AC2_BD2_dis])

    if max_side == AC1_BD1_dis:
        A = unordered[0]
        B = unordered[2]
        C = unordered[1]
        D = unordered[3]
    if max_side == AC1_BD2_dis:
        A = unordered[0]
        B = unordered[3]
        C = unordered[1]
        D = unordered[2]
    if max_side == AC2_BD1_dis:
        A = unordered[1]
        B = unordered[2]
        C = unordered[0]
        D = unordered[3]
    if max_side == AC2_BD2_dis:
        A = unordered[1]
        B = unordered[3]
        C = unordered[0]
        D = unordered[2]

    vertices_out = np.array([A,B,C,D])
    return vertices_out

def embedded_vert(vertices):
    """
    Return the vertices embedded in the original vertices

    The input is four coordinates of a quadrangle in order ABCD. The output is
    the four coordinates on the centres of the sides like: [AB,BC,CD,DA]

    Parameters
    ----------
    vertices : numpy array of 4 x 3
        Four coordinates [A,B,C,D]
    Returns
    -------
    out_vertices : numpy array of 4 x 3
        Four embedded coordinates [AB,BC,CD,DA]

    """
    out_vertices = []

    out_vertices.append(np.mean([vertices[0],vertices[1]],axis=0))
    out_vertices.append(np.mean([vertices[1],vertices[2]],axis=0))
    out_vertices.append(np.mean([vertices[2],vertices[3]],axis=0))
    out_vertices.append(np.mean([vertices[3],vertices[0]],axis=0))

    out_vertices = np.array(out_vertices)

    return out_vertices

def project_point(point, plane_coeffs):
    """
    Project a point onto a plane

    Parameters
    ----------
    point : length 3 numpy array
        The point to be projected on the plane
    plane_coeffs : length 4 numpy array
        Plane equation coefficients such that a point on the plane is:
        ax + by + cz + d = 0. The array is [a,b,c,d]
    Returns
    -------
    out_point : length 3 numpy array
        The coordinates of the projected point

    """
    plane_normal_vector = plane_coeffs[:3]
    point_norm2 = np.sum(plane_normal_vector*plane_normal_vector)
    proj_parameter = plane_normal_vector * point / point_norm2 + plane_coeffs[3]

    projection = point - plane_normal_vector * proj_parameter

    return projection

def project_pair_to_vector(coord_pair, plane_coeffs):
    """
    Return the normalised vector formed by a coordinate pair projected onto a plane

    Parameters
    ----------
    coord_pair : numpy array of 2 x 3
        Pair of coordinates
    plane_coeffs : length r numpy array
        Plane equation coefficients such that a point on the plane is:
        ax + by + cz + d = 0. The array is [a,b,c,d]
    Returns
    -------
    vector : length 3 numpy array
        The normalised vector resulting from the projection of the pair of
        points

    """
    proj_a = project_point(coord_pair[0], plane_coeffs)
    proj_b = project_point(coord_pair[1], plane_coeffs)

    vector = proj_b - proj_a

    vector /= np.linalg.norm(vector)
    return vector

def project_quad_to_vectors(quad_vert, plane_coeffs):
    """
    Return the normalised vectors formed by the diagonals of a quadrangle

    Parameters
    ----------
    quad_vert : numpy array 4 x 3
        Four coordinates [A,B,C,D]
    plane_coeffs : length 4 numpy array
        Plane equation coefficients such that a point on the plane is:
        ax + by + cz + d = 0. The array is [a,b,c,d]
    Returns
    -------
    vectors : length 2 x 3 numpy array
        The normalised vectors resulting from the projection of the pairs of
        points A -> C and B -> D

    """
    vec_1 = project_pair_to_vector([quad_vert[0],quad_vert[2]],plane_coeffs)
    vec_2 = project_pair_to_vector([quad_vert[1],quad_vert[3]],plane_coeffs)

    vectors = np.array([vec_1,vec_2])

    return vectors

def vec_angle(vector_1, vector_2, degrees = True):
    """
    Return the angle between two numpy vectors.

    The angle is brought into the range [-180,180] for degrees or [-1,1] for
    radians

    Parameters
    ----------
    vector_1, vector_2 : N x 1 numpy array
        The vectors whose angle needs to be calculated.
    degrees : bool (optional)
        Result in degrees or in radians. Default = True, so degrees
    Returns
    -------
    out_angle : float
        The angle between vectors
    """
    norm_1 = np.linalg.norm(vector_1)
    norm_2 = np.linalg.norm(vector_2)
    dot = np.dot(vector_1,vector_2)

    ang = np.arccos(dot/(norm_1*norm_2)) % (2 * np.pi)

    if degrees:
        ang = np.degrees(ang)
    return ang
