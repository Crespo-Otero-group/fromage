import numpy as np
from pytest import approx

def test_complete_mol(hc1_cell):
    new_mol, new_cell = hc1_cell.complete_mol(0)
    sel = new_mol.select(0)
    assert len(new_mol) == len(sel)

def test_complete_cell(hc1_cell):
    new_cell, new_mols = hc1_cell.complete_cell()
    assert len(new_mols[0]) == 37

def test_supercell(hc1_cell):
    trans = np.array([2,2,2])
    new_cell = hc1_cell.supercell(trans)
    assert len(new_cell) == 1184

def test_big_supercell(hc1_cell):
    trans = np.array([3,3,3])
    new_cell = hc1_cell.supercell(trans)
    assert len(new_cell) == 3996

def test_centered_supercell(hc1_cell):
    trans = np.array([1,1,1])
    new_cell = hc1_cell.centered_supercell(trans)
    assert len(new_cell) == approx(3996)

def test_centered_supercell_alt(hc1_cell):
    trans = np.array([1,1,1])
    new_cell = hc1_cell.centered_supercell(trans, from_origin=True)
    assert len(new_cell) == approx(1184)

def test_make_cluster(hc1_cell):
    clust = hc1_cell.make_cluster(10)
    assert len(clust) == 74

def test_confine(hc1_complete_cell):
    conf = hc1_complete_cell.confined()
    assert conf[19].x == approx(-0.202155)
