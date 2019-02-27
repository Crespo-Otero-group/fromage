import numpy as np


def test_coord_array(h2o_dimer):
    arr = h2o_dimer.coord_array()
    print(arr)
    from sklearn.metrics.pairwise import pairwise_distances
    print(pairwise_distances(arr))
    assert np.shape(arr) == (6, 3)

def test_pairwise_distances2(h2o_dimer):
    dis = h2o_dimer.pairwise_distances2()
