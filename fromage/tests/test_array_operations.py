import fromage.utils.array_operations as ao
import numpy as np
from pytest import approx

def test_dist_mat(h2o_dim_dist_arr):
    dis = ao.dist_mat(h2o_dim_dist_arr)
    assert dis[1][0] == approx(0.91,rel=0.1)


def test_find_largest(h2o_dim_dist_arr):
    largest_inds = ao.find_largest(h2o_dim_dist_arr,3)
    assert largest_inds == [(4, 0), (4, 2), (5, 0)]

def test_plane_from_coord(hc1_array):
    results = ao.plane_from_coord(hc1_array)
    arr = np.array(results)
    desired = np.array([-0.29188565799088845, 0.742637122670864, -0.6027378092424318, 1.7825747834180873e-07])
    np.testing.assert_allclose(arr,desired)

def test_extreme_pairs(h2o_dim_array):
    res = ao.extreme_pairs(h2o_dim_array,2)
    check = np.array([[[ 3.758602,0.5,0.504284],[ 0.,0.,0.,]],
        [[ 3.758602,0.5,0.504284],[ 0.260455,0.,-0.872893]]])
    assert (res == check).all()
