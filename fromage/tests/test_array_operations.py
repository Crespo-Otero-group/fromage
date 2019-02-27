import fromage.utils.array_operations as ao
def test_find_largest(h2o_dim_dist_arr):
    largest_inds = ao.find_largest(h2o_dim_dist_arr,3)
    assert largest_inds == [(4, 0), (4, 2), (5, 0)]
    return

