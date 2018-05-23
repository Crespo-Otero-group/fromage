import pytest
import cryspy.io.read_file as rf
from cryspy.utils.atom import Atom
from cryspy.utils.mol import Mol

@pytest.fixture
def at_list():
    """Return a water dimer atom list"""
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
def hc1_quad():
    """Return a HC1 quadrimer"""
    out_mo = Mol(rf.read_xyz("hc1_quad.xyz")[-1])
    return out_mo

@pytest.fixture
def newat():
    """Return an Atom object C at origin"""
    return Atom("C", 0.0, 0.0, 0.0)

@pytest.fixture
def empty_mol():
    return Mol([])

def test_len(h2o_dimer):
    """The len method is implemented"""
    assert len(h2o_dimer) == 6

def test_len_empty(empty_mol):
    assert len(empty_mol) == 0

def test_for_loop(h2o_dimer):
    counter = 0
    for at in h2o_dimer:
        counter += 1
    assert counter == 6

def test_add(h2o_dimer, newat):
    """Addition is implemented"""
    new_mol = h2o_dimer + Mol(newat)
    assert len(new_mol) == 7

def test_append(h2o_dimer, newat):
    """Appending is implemented"""
    h2o_dimer.append(newat)
    assert len(h2o_dimer) == 7

def test_select_h2o_dimer(h2o_dimer):
    water = h2o_dimer.select(3)
    assert len(water) == 3

def test_select_hc1_quad(hc1_quad):
    mol = hc1_quad.select([0])
    assert len(mol) == 37

def test_multiselect_hc1_quad(hc1_quad):
    mol = hc1_quad.select([0,74])
    mol.write_xyz("hi.xyz")
    assert len(mol) == 74

