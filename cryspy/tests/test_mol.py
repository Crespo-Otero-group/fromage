import pytest
import cryspy.io.read_file as rf
from cryspy.utils.atom import Atom
from cryspy.utils.mol import Mol

@pytest.fixture
def at_list():
    """Returns a water dimer atom list"""
    out_list = rf.read_pos("h2o_dimer.xyz")
    return out_list

def test_at_list_type(at_list):
    """The atom list is made of Atom objects"""
    assert all(isinstance(i, Atom) for i in at_list)

def test_init():
    """The Mol object is initialised correctly"""
    mo = Mol(at_list())
    assert isinstance(mo,Mol)

@pytest.fixture
def h2o_dimer(at_list):
    """Return a water dimer Mol object"""
    out_mo = Mol(at_list)
    return out_mo

@pytest.fixture
def newat():
    """Return an Atom object C at origin"""
    return Atom("C", 0.0, 0.0, 0.0)

def test_len(h2o_dimer):
    """The len method is implemented"""
    assert len(h2o_dimer) == 6

def test_add(h2o_dimer, newat):
    """Addition is implemented"""
    new_mol = h2o_dimer + Mol(newat)
    assert len(new_mol) == 7

def test_append(h2o_dimer, newat):
    """Appending is implemented"""
    h2o_dimer.append(newat)
    assert len(h2o_dimer) == 7

