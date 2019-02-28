import numpy as np
from pytest import approx

def test_coord_array(h2o_dimer):
    arr = h2o_dimer.coord_array()
    assert np.shape(arr) == (6, 3)

def test_extreme_at_pairs(h2o_dimer):
    res = h2o_dimer.extreme_at_pairs(2)
    check = np.array([[[ 3.758602,0.5,0.504284],[ 0.,0.,0.,]],
        [[ 3.758602,0.5,0.504284],[ 0.260455,0.,-0.872893]]])
    assert (res == check).all()
