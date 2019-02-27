import fromage.utils.array_operations as ao
from pytest import approx

def test_dist_mat(h2o_dim_dist_arr):
    dis = ao.dist_mat(h2o_dim_dist_arr)
    assert dis[1][0] == approx(0.91,rel=0.1)


def test_find_largest(h2o_dim_dist_arr):
    largest_inds = ao.find_largest(h2o_dim_dist_arr,3)
    assert largest_inds == [(4, 0), (4, 2), (5, 0)]

def test_plane_from_coord(hc1_array):
    print(ao.plane_from_coord(hc1_array))
