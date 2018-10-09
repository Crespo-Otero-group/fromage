import pytest
from pytest import approx
import numpy as np

from cryspy.utils.atom import Atom


@pytest.fixture
def c_at():
    """Return a C Atom object at origin"""
    out_at = Atom("C", 0.0, 0.0, 0.0)
    return out_at


@pytest.fixture
def o_at():
    """Return am O Atom object at x = 0.8"""
    out_at = Atom("O", 0.8, 0.0, 0.0)
    return out_at

@pytest.fixture
def o_at_outside():
    """Return am O Atom object outside of the 111 cell"""
    out_at = Atom("O", 2.1, -0.3, 1.4)
    return out_at

@pytest.fixture
def vectors():
    vec_out = np.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]])
    return vec_out


def test_dist_to_atom(c_at, o_at):
    """Distance to atom is implemented"""
    assert c_at.dist(o_at) == approx(0.8)


def test_at_lap(c_at, o_at):
    """Atom vdw overlap works"""
    assert c_at.dist(o_at,ref='vdw') == approx(-2.42)


def test_per_dist(c_at, o_at, vectors):
    assert c_at.per_dist(o_at, vectors) == approx(0.2)


def test_per_lap(c_at, o_at, vectors):
    assert c_at.per_dist(o_at, vectors, ref='vdw') == approx(-3.02)


def test_per_lap_img(c_at, o_at, vectors):
    r, new_at = c_at.per_dist(o_at, vectors, ref='vdw', new_pos=True)
    assert o_at != new_at

def test_put_in_cell(o_at_outside,vectors):
    in_box = o_at_outside.put_in_cell(vectors)
    assert in_box.x == approx(0.1)
