import pytest
import fromage.io.read_file as rf
import numpy as np
from fromage.utils.mol import Mol
from fromage.utils.atom import Atom

# Atom fixtures
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


# Mol fixtures
@pytest.fixture
def at_list():
    """Return a water dimer atom list"""
    out_list = rf.read_pos("h2o_dimer.xyz")
    return out_list

@pytest.fixture
def h2o_dimer(at_list):
    """Return a water dimer Mol object"""
    out_mo = Mol(at_list)
    return out_mo


@pytest.fixture
def hc1_quad():
    """Return a HC1 quadrimer"""
    out_mo = Mol(rf.read_pos("hc1_quad.xyz"))
    return out_mo


@pytest.fixture
def hc1_cell():
    """Return an HC1 cell"""
    cell = Mol(rf.read_pos("hc1_cell.xyz"))
    vectors = np.array([[12.1199998856, 0.0, 0.0],
                        [0.0, 10.2849998474, 0.0],
                        [-5.4720203118, 0.0, 11.2441994632]])
    cell.vectors = vectors
    return cell

@pytest.fixture
def hc1_complete_cell():
    """Return an HC1 cell"""
    cell = Mol(rf.read_pos("hc1_complete_cell.xyz"))
    vectors = np.array([[12.1199998856, 0.0, 0.0],
                        [0.0, 10.2849998474, 0.0],
                        [-5.4720203118, 0.0, 11.2441994632]])
    cell.vectors = vectors
    return cell

@pytest.fixture
def newat():
    """Return an Atom object C at origin"""
    return Atom("C", 0.0, 0.0, 0.0)

@pytest.fixture
def empty_mol():
    return Mol([])

@pytest.fixture
def c_o():
    c_at = Atom("C", 0.0, 0.0, 0.0)
    o_at = Atom("O", 1.0, 0.0, 0.0)
    out_mol = Mol([c_at, o_at])
    return out_mol

@pytest.fixture
def h2o_dup():
    out_mol = rf.mol_from_file("h2o_repeated.xyz")
    return out_mol


