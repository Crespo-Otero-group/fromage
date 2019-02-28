import numpy as np
from pytest import approx
from numpy.testing import assert_allclose

def test_coord_array(h2o_dimer):
    arr = h2o_dimer.coord_array()
    assert np.shape(arr) == (6, 3)

def test_calc_coord_array(h2o_dimer):
    h2o_dimer.calc_coord_array()
    assert np.shape(h2o_dimer.geom_info.coord_array) == (6,3)

def test_plane_coeffs(hc1_mol):
    coeffs = hc1_mol.plane_coeffs()
    expected = np.array([-2.91885658e-01,7.42637123e-01,-6.02737809e-01,1.78257478e-07])
    assert_allclose(coeffs, expected)


#def test_calc_axes(hc1_mol):
#    hc1_mol.calc_axes()
#    print(hc1.prin_ax)
#    print(hc1.sec_ax)
