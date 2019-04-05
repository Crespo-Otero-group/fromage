"""The functions having to do with planes and quadrangles projected thereupon"""
import numpy as np

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
    from fromage.utils.array_operations import dist_mat, find_largest, distance
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

    t = (-plane_coeffs[3]-plane_coeffs[0]*point[0]-plane_coeffs[1]*point[1]-plane_coeffs[2]*point[2])
    projection = np.array([point[0]+t*plane_coeffs[0],point[1]+t*plane_coeffs[1],point[2]+t*plane_coeffs[2]])

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
    Return the normalised vectors, projected onto a plane, formed by the
    diagonals of a quadrangle

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
