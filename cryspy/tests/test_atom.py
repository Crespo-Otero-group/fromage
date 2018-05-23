import pytest
import cryspy.io.read_file as rf
from pytest import approx
from cryspy.utils.atom import Atom
from cryspy.utils.mol import Mol

@pytest.fixture
def c_at():
    """Return a C Atom object at origin"""
    out_at = Atom("C", 0.0, 0.0, 0.0)
    return out_at

@pytest.fixture
def h_at():
    """Return am O Atom object at x = 1.3"""
    out_at = Atom("O", 1.3, 0.0, 0.0)
    return out_at

def test_dist_to_atom(c_at, h_at):
    """Distance to atom is implemented"""
    assert c_at.dist_at(h_at) == approx(1.3)

def test_at_lap(c_at, h_at):
    """Atom vdw overlap works"""
    assert c_at.at_lap(h_at) == approx(1.92)
