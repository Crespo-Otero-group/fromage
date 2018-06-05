import pytest
from pytest import approx
import cryspy.io.read_file as rf
import numpy as np
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
    assert isinstance(mo, Mol)


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
    mol = hc1_quad.select([0, 74])
    assert len(mol) == 74


def test_per_select_hc1_cell(hc1_cell):
    """Check that periodic select gets all atoms"""
    selected = hc1_cell.per_select(0)
    assert len(selected) == 37


def test_per_select_complete(hc1_cell, hc1_quad):
    """Check that periodic select completes the molecule"""
    selected = hc1_cell.per_select(0)
    new, old = hc1_cell.per_select(0, old_pos=True)
    new_sel = selected.select(0)
    assert len(selected) == len(new_sel)


def test_multi_per_select(hc1_cell, hc1_quad):
    selected = hc1_cell.per_select([0, 1])
    new_sel = hc1_quad.select(0)
    assert len(selected) == len(new_sel) * 2


def test_mol_segregation(hc1_quad):
    mols = hc1_quad.segregate()
    assert len(mols) == 4 and len(mols[1]) == 37


def test_complete_mol(hc1_cell):
    new_mol, new_cell = hc1_cell.complete_mol(0)
    sel = new_mol.select(0)
    assert len(new_mol) == len(sel)

def test_complete_cell(hc1_cell):
    new_cell, new_mols = hc1_cell.complete_cell()
    new_cell.write_xyz("hi.xyz")
    assert len(new_mols[0]) == 37

def test_centroid(c_o):
    cen = c_o.centroid()
    assert cen[0] == approx(0.5)

def test_center(c_o):
    c_o.center_mol()
    assert c_o[0].x == approx(-0.5)

def test_translate(c_o):
    vec = np.array([0.25,1.0,0.0])
    c_o.translate(vec)
    assert c_o[1].x == approx(1.25)

def test_supercell(hc1_cell):
    trans = np.array([2,2,2])
    new_cell = hc1_cell.supercell(trans)
    assert len(new_cell) == 1184

def test_centered_supercell(hc1_cell):
    trans = np.array([1,1,1])
    new_cell = hc1_cell.centered(trans)
    assert len(new_cell) == approx(1332)

def test_centered_supercell(hc1_cell):
    trans = np.array([1,1,1])
    new_cell = hc1_cell.centered(trans, from_origin=True)
    assert len(new_cell) == approx(592)
