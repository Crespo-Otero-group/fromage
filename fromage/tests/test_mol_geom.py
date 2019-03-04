import numpy as np
from numpy.testing import assert_allclose

def test_coord_array(h2o_dimer):
    arr = h2o_dimer.coord_array()
    assert np.shape(arr) == (6, 3)

def test_calc_coord_array(h2o_dimer):
    h2o_dimer.calc_coord_array()
    assert np.shape(h2o_dimer.geom.coord_array) == (6,3)

def test_calc_plane_coeffs(hc1_mol):
    hc1_mol.calc_plane_coeffs()
    coeffs = hc1_mol.geom.plane_coeffs
    expected = np.array([-2.91885658e-01,7.42637123e-01,-6.02737809e-01,1.78257478e-07])
    assert_allclose(coeffs, expected)

def test_calc_axes_hc1(hc1_mol):
    hc1_mol.calc_axes()
    expected_princ = np.array([0.970348,0.062753,-0.233425])
    expected_sec = np.array([0.183283,0.438568,0.87981])
    assert_allclose(hc1_mol.geom.prin_ax,expected_princ, rtol=1e-4)
    assert_allclose(hc1_mol.geom.sec_ax,expected_sec, rtol=1e-4)

def test_calc_axes(rectangle_mol):
    rectangle_mol.calc_axes()
    expected_princ = np.array([0.999931,0.011769,0.])
    expected_sec = np.array([0.011769,-0.999931,0.])
    expected_perp = np.array([ 0., 0.,  -1.      ])
    assert_allclose(rectangle_mol.geom.prin_ax,expected_princ, rtol=1e-4)
    assert_allclose(rectangle_mol.geom.sec_ax,expected_sec, rtol=1e-4)
    assert_allclose(rectangle_mol.geom.perp_ax,expected_perp, rtol=1e-4)
