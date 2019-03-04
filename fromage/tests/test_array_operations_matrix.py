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
