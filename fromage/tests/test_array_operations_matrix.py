import fromage.utils.array_operations as ao
import numpy as np
from numpy.testing import assert_allclose
from pytest import approx

def test_cross_product_matrix(vec_111):
    res = ao.cross_product_matrix(vec_111)
    expected = np.array([[0.,-1.,1.],
                         [1.,0.,-1.],
                         [-1.,1.,0.]])
    assert_allclose(res,expected)

def test_rotation_matrix(vec_100):
    mat = ao.rotation_matrix(vec_100,90)
    res = np.dot(mat,np.array([0.,1.,0.]))
    expected = np.array([0.,0.,1.])
    assert_allclose(res,expected,atol=1e-7)

def test_rotation_matrix_radians(vec_100):
    mat = ao.rotation_matrix(vec_100,np.pi/2,degrees=False)
    res = np.dot(mat,np.array([0.,1.,0.]))
    expected = np.array([0.,0.,1.])
    assert_allclose(res,expected,atol=1e-7)
