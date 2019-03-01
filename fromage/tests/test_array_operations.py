import fromage.utils.array_operations as ao
import numpy as np
from numpy.testing import assert_allclose
from pytest import approx

def test_distance(vec_100,vec_220):
    distance = ao.distance(vec_100,vec_220)
    assert distance == approx(2.236,rel=0.001)

def test_closest(vec_100,vec_220):
    origin = np.array([0.,0.,0.])
    close_point, dis = ao.closest(origin,[vec_100,vec_220])
    assert_allclose(close_point,vec_100)
    assert dis == approx(1.)

def test_dist_mat(h2o_dim_dist_arr):
    dis = ao.dist_mat(h2o_dim_dist_arr)
    assert dis[1][0] == approx(0.91,rel=0.1)

def test_find_largest(h2o_dim_dist_arr):
    largest_inds = ao.find_largest(h2o_dim_dist_arr,3)
    assert largest_inds == [(4, 0), (4, 2), (5, 0)]

