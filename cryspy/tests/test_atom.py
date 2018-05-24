import pytest
from pytest import approx
import numpy as np

import cryspy.io.read_file as rf
from cryspy.utils.atom import Atom
from cryspy.utils.mol import Mol

@pytest.fixture
def c_at():
    """Return a C Atom object at origin"""
    out_at = Atom("C", 0.0, 0.0, 0.0)
    return out_at

@pytest.fixture
def o_at():
    """Return am O Atom object at x = 1.3"""
    out_at = Atom("O", 1.3, 0.0, 0.0)
    return out_at

@pytest.fixture
def vectors():
    vec_out = np.array([[1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])
    return vec_out

def test_dist_to_atom(c_at, o_at):
    """Distance to atom is implemented"""
    assert c_at.dist_at(o_at) == approx(1.3)

def test_at_lap(c_at, o_at):
    """Atom vdw overlap works"""
    assert c_at.at_lap(o_at) == approx(1.92)

def test_per_dist(c_at, o_at, vectors):
    assert c_at.per_dist(o_at,vectors) == approx(0.3)

def test_per_lap(c_at, o_at, vectors):
    assert c_at.per_lap(o_at,vectors) == approx(2.92)

def test_per_lap_img(c_at, o_at, vectors):
    r, new_at = c_at.per_lap(o_at,vectors,new_pos=True)
    assert o_at != new_at
