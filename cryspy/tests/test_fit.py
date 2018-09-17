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

@pytest.fixture
def benz_pot_cub():
    """Return the CubeGrid for a benzene cell potential"""
    cub, mol = rf.read_cube("benz_pot.cube")
    return cub

def test_fit_benz_clust(benz_clust_char):
    sample_point = Atom("",0,0,0)
    sample_point2 = Atom("",0.1,0,0)
    sample_point.es = 0
    sample_point2.es = 0
    fi.fit_points(benz_clust_char, None, Mol([sample_point, sample_point2]))
    return

def test_fit_benz_clust_shell(benz_clust_char):
    sample_point = Atom("",0,0,0)
    sample_point2 = Atom("",0.1,0,0)
    sample_point.es = 0
    sample_point2.es = 0
    fixed_atoms = benz_clust_char.select(0)
    var_atoms = Mol([i for i in benz_clust_char if i not in fixed_atoms])
    fi.fit_points(var_atoms, fixed_atoms,Mol([sample_point, sample_point2]))
    return

def test_fit_to_pot(benz_pot_cub,benz_clust_char):
    print(benz_pot_cub.grid)
    sample_point = Atom("",0,0,0)
    sample_point.es = 0.17671
    print(benz_clust_char.es_pot([0,0,0]))
    fixed_atoms = benz_clust_char.select(0)
    var_atoms = Mol([i for i in benz_clust_char if i not in fixed_atoms])
    fi.fit_points(var_atoms, fixed_atoms,Mol([sample_point]))
    out_clust = var_atoms + fixed_atoms
    print(out_clust.es_pot([0,0,0]))
    return

def test_sample(benz_pot_cub,benz_cell):
    pts = benz_pot_cub.shell_region(benz_cell, 0.2, 0.7)
    out_mol = Mol([])
    for point in pts:
        out_mol.append(Atom("point",point[0],point[1],point[2]))
    out_mol.write_xyz("tmp.xyz")
