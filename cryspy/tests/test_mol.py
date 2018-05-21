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
    assert all(isinstance(i, Atom) for i in at_list)

def test_init():
    mo = Mol(at_list())
    assert isinstance(mo,Mol)

@pytest.fixture
def mo(at_list):
    """Returns a water dimer Mol object"""
    out_mo = Mol(at_list)
    return out_mo

def test_iter(mo):
    for at in mo:
        print(at)
    print(mo)

