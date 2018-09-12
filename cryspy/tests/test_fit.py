import pytest
from pytest import approx
import cryspy.io.read_file as rf
import cryspy.utils.fit as fi
import cryspy.scripts.assign_charges as ac
from cryspy.utils.atom import Atom
from cryspy.utils.mol import Mol
import numpy as np
from copy import deepcopy

@pytest.fixture
def benz_cell():
    """Return a Mol object of the uncharged benzene cell"""
    out_cell = rf.mol_from_file("benzene_cell.xyz")
    out_cell.vectors = rf.read_vectors("benzene_vectors")
    return out_cell

@pytest.fixture
def benz_solo():
    """Return a Mol object of a charged benzene"""
    out_mol = rf.mol_from_gauss("benz.log")
    return out_mol

@pytest.fixture
def benz_cell_char(benz_cell,benz_solo):
    """Return a Mol object of a charged benzene cell"""
    out_char_cell = deepcopy(benz_cell)
    ac.assign_charges(benz_solo, None, out_char_cell, out_char_cell.vectors, 1.7)
    return out_char_cell

@pytest.fixture
def benz_clust_char(benz_solo):
    """Return a Mol object of a charged benzene cluster"""
    out_char_clust = rf.mol_from_file("benz_clust.xyz")
    ac.assign_charges(benz_solo, None, out_char_clust, None, 1.7)
    return out_char_clust

def test_fit_benz_clust(benz_clust_char):
    sample_point = Atom("",0,0,0)
    sample_point2 = Atom("",0.1,0,0)
    sample_point.es = 0
    sample_point2.es = 0
    fi.fit_points(benz_clust_char, Mol([sample_point, sample_point2]))

    return
