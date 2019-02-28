import pytest
import fromage.io.read_file as rf
import numpy as np
import fromage.utils.array_operations as ao
from fromage.utils.mol import Mol
from fromage.utils.atom import Atom
from scipy.spatial.distance import cdist

# Atom fixtures
@pytest.fixture
def c_at():
    """C Atom object at origin"""
    out_at = Atom("C", 0.0, 0.0, 0.0)
    return out_at


@pytest.fixture
def o_at():
    """O Atom object at x = 0.8"""
    out_at = Atom("O", 0.8, 0.0, 0.0)
    return out_at

@pytest.fixture
def o_at_outside():
    """O Atom object outside of the 111 cell"""
    out_at = Atom("O", 2.1, -0.3, 1.4)
    return out_at

@pytest.fixture
def vectors():
    """Eye 3x3 vectors"""
    vec_out = np.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]])
    return vec_out


# Mol fixtures
@pytest.fixture
def at_list():
    """Water dimer atom list"""
    out_list = rf.read_pos("h2o_dimer.xyz")
    return out_list

@pytest.fixture
def h2o_dimer(at_list):
    """Water dimer Mol object"""
    out_mo = Mol(at_list)
    return out_mo


@pytest.fixture
def hc1_quad():
    """HC1 quadrimer"""
    out_mo = Mol(rf.read_pos("hc1_quad.xyz"))
    return out_mo


@pytest.fixture
def hc1_cell():
    """HC1 cell"""
    cell = Mol(rf.read_pos("hc1_cell.xyz"))
    vectors = np.array([[12.1199998856, 0.0, 0.0],
                        [0.0, 10.2849998474, 0.0],
                        [-5.4720203118, 0.0, 11.2441994632]])
    cell.vectors = vectors
    return cell

@pytest.fixture
def hc1_complete_cell():
    """HC1 completed cell"""
    cell = Mol(rf.read_pos("hc1_complete_cell.xyz"))
    vectors = np.array([[12.1199998856, 0.0, 0.0],
                        [0.0, 10.2849998474, 0.0],
                        [-5.4720203118, 0.0, 11.2441994632]])
    cell.vectors = vectors
    return cell

@pytest.fixture
def newat():
    """Atom object C at origin"""
    return Atom("C", 0.0, 0.0, 0.0)

@pytest.fixture
def empty_mol():
    """Empty Mol object"""
    return Mol([])

@pytest.fixture
def c_o():
    """CO Mol object"""
    c_at = Atom("C", 0.0, 0.0, 0.0)
    o_at = Atom("O", 1.0, 0.0, 0.0)
    out_mol = Mol([c_at, o_at])
    return out_mol

@pytest.fixture
def h2o_dup():
    """H2O molecule with duplicate atoms"""
    out_mol = rf.mol_from_file("h2o_repeated.xyz")
    return out_mol

# Array fixtures
@pytest.fixture
def h2o_dim_array():
    """Coordinate array for H2O dimer"""
    mol = rf.mol_from_file("h2o_dimer.xyz")
    arr = mol.coord_array()
    return arr

@pytest.fixture
def h2o_dim_dist_arr(h2o_dim_array):
    """Distance array for the H2O dimer"""
    return(ao.dist_mat(h2o_dim_array))

@pytest.fixture
def hc1_array():
    """Coordinate array for the HC1 monomer"""
    mol = rf.mol_from_file("hc1_mol.xyz")
    arr = mol.coord_array()
    return arr

@pytest.fixture
def rectangle_pairs_array():
    """2 pairs of coordinates defining the 4 X 2 rectangle"""
    lis_rectangle = [[[0.,0.,0.],[4.,2.,0.]],[[0.,2.,0.],[4.,0.,0.]]]
    arr = np.array(lis_rectangle)
    return arr
