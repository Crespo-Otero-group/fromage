import fromage.utils.array_operations as ao
import numpy as np
from numpy.testing import assert_allclose
from pytest import approx

# Angle tests
def test_angle_between_vectors(vec_100, vec_220):
    ang = ao.vec_angle(vec_100, vec_220)
    assert ang == approx(45.)

def test_angle_between_vectors_rad(vec_100, vec_220):
    ang = ao.vec_angle(vec_100, vec_220, degrees = False)
    assert ang == approx(45.*np.pi/180)
