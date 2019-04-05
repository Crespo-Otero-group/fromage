import fromage.utils.array_operations as ao
import numpy as np
from numpy.testing import assert_allclose
from pytest import approx

def test_plane_from_coord(rectangle_array):
    results = ao.plane_from_coord(rectangle_array)
    arr = results
    expected = np.array([0.,0.,1.,0.])
    assert_allclose(arr,expected)

def test_plane_from_coord2(hc1_array):
    results = ao.plane_from_coord(hc1_array)
    arr = results
    expected = np.array([-0.29188565799088845, 0.742637122670864, -0.6027378092424318, 1.7825747834180873e-07])
    assert_allclose(arr,expected)

# Tweaked flat rectangle tests
def test_arb_quad_from_coord(arbitrary_flat_points, arbitrary_flat_vertices):
    res = ao.quadrangle_from_coord(arbitrary_flat_points)
    assert_allclose(res,arbitrary_flat_vertices)

def test_arb_embedded_vert(arbitrary_flat_vertices):
    new_pairs = ao.embedded_vert(arbitrary_flat_vertices)
    expected = np.array([[3.  , 3.05, 0.  ],
                         [6.1  , 1.45, 0.  ],
                         [3.1, -0.1, 0.],
                         [0. , 1.5, 0. ]])
    assert_allclose(new_pairs,expected)

def test_arb_project_quad_to_vectors(arbitrary_flat_vertices,z_plane_coeffs):
    emb_vert = ao.embedded_vert(arbitrary_flat_vertices)
    projected_vecs = ao.project_quad_to_vectors(emb_vert,z_plane_coeffs)
    expected = np.array([[ 0.03173 , -0.999496,  0.      ],
                        [-0.999966,  0.008196,  0.      ]])
    assert_allclose(projected_vecs, expected, rtol=1e-4)
    return

# Plain rectangle tests
def test_embedded_vert(rectangle_array):
    new_vert = ao.embedded_vert(rectangle_array)
    expected = np.array([[0.,1.,0.],[2.,2.,0.],[4.,1.,0.],[2.,0.,0.]])
    assert_allclose(new_vert,expected)

def test_project_on_plane(z_plane_coeffs):
    point = np.array([1.,1.,2.5])
    projected = ao.project_point(point,z_plane_coeffs)
    expected = np.array([1.,1.,0.])
    assert_allclose(projected,expected)

# Random rectangle tests
def test_project_pair_to_vector(arbitrary_pair,z_plane_coeffs):
    projected_vec = ao.project_pair_to_vector(arbitrary_pair,z_plane_coeffs)
    expected = np.array([1.,0.,0.])
    assert_allclose(projected_vec, expected)
    return

def test_project_quad_to_vectors(arbitrary_vertices,z_plane_coeffs):
    projected_vecs = ao.project_quad_to_vectors(arbitrary_vertices,z_plane_coeffs)
    expected = np.array([[1.,0.,0.],[0.,-1.,0.]])
    assert_allclose(projected_vecs, expected)
    return

# Triangle set tests
def test_tri_quad_from_coord(triangle_shape_coord,triangle_corners_4):
    res = ao.quadrangle_from_coord(triangle_shape_coord)
    assert_allclose(res,triangle_corners_4)

def test_tri_embedded_vert(triangle_corners_4):
    new_pairs = ao.embedded_vert(triangle_corners_4)
    expected = np.array([[3.5, 3., 0.],
                         [6.,2.,0.],
                         [3.5, 1.5, 0.],
                         [1.,2.5,0.]])
    assert_allclose(new_pairs,expected)

def test_tri_project_quad_to_vectors(triangle_corners_4,z_plane_coeffs):
    emb_vert = ao.embedded_vert(triangle_corners_4)
    projected_vecs = ao.project_quad_to_vectors(emb_vert,z_plane_coeffs)
    expected = np.array([[ 0., -1.,  0.      ],
                        [-0.995037,  0.099504,  0.      ]])
    assert_allclose(projected_vecs, expected, rtol=1e-4)
    return

# Angle tests
def test_angle_between_vectors(vec_100, vec_220):
    #print(ao.vec_angle(vec_100, vec_220))
    return
