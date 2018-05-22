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
    assert c_at.dist_at(h_at) == 1.3

def test_smart_atom_distance(c_at, h_at):
    """Smart distance is implemented"""
    approx(c_at.smart_dist_at(h_at), 1.92)
