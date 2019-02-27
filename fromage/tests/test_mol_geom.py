import numpy as np
from pytest import approx

def test_coord_array(h2o_dimer):
    arr = h2o_dimer.coord_array()
    assert np.shape(arr) == (6, 3)

def test_pairwise_distances(h2o_dimer):
    dis = h2o_dimer.pairwise_distances()
    print(dis)
    assert dis[0][1] == approx(0.91,rel=0.1)
