import fromage.utils.array_operations as ao
import numpy as np
from numpy.testing import assert_allclose
from pytest import approx

def test_dist_mat(h2o_dim_dist_arr):
    dis = ao.dist_mat(h2o_dim_dist_arr)
    assert dis[1][0] == approx(0.91,rel=0.1)

def test_find_largest(h2o_dim_dist_arr):
    largest_inds = ao.find_largest(h2o_dim_dist_arr,3)
    assert largest_inds == [(4, 0), (4, 2), (5, 0)]

def test_plane_from_coord(rectangle_array):
    results = ao.plane_from_coord(rectangle_array)
    arr = results
    expected = np.array([0.,0.,1.,0.])
    assert_allclose(arr,expected)

def test_plane_from_coord(hc1_array):
    results = ao.plane_from_coord(hc1_array)
    arr = results
    expected = np.array([-0.29188565799088845, 0.742637122670864, -0.6027378092424318, 1.7825747834180873e-07])
    assert_allclose(arr,expected)

def test_extreme_pairs(h2o_dim_array):
    res = ao.extreme_pairs(h2o_dim_array,2)
    expected = np.array([[[ 3.758602,0.5,0.504284],[ 0.,0.,0.,]],
        [[ 3.758602,0.5,0.504284],[ 0.260455,0.,-0.872893]]])
    assert_allclose(res,expected)

def test_embedded_pairs(rectangle_pairs_array):
    new_pairs = ao.embedded_pairs(rectangle_pairs_array)
    expected = np.array([[[0.,1.,0.],[4.,1.,0.]],[[2.,0.,0.],[2.,2.,0.]]])
    assert_allclose(new_pairs,expected)

def test_project_on_plane(rectangle_array,rectangle_plane_coeffs):
    point = np.array([1.,1.,2.5])
    projected = ao.project_point(point,rectangle_plane_coeffs)
    expected = np.array([1.,1.,0.])
    assert_allclose(projected,expected)

def test_project_pair_to_vector(arbitrary_vector,rectangle_plane_coeffs):
    projected_vec = ao.project_pair_to_vector(arbitrary_vector,rectangle_plane_coeffs)
    expected = np.array([-1.,0.,0.])
    assert_allclose(projected_vec, expected)
    return

def test_project_pairs_to_vectors(arbitrary_pairs_array,rectangle_plane_coeffs):
    projected_vecs = ao.project_pairs_to_vectors(arbitrary_pairs_array,rectangle_plane_coeffs)
    expected = np.array([[-1.,0.,0.],[0.,1.,0.]])
    assert_allclose(projected_vecs, expected)
    return
