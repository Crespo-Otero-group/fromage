import os
import pytest
import fromage.io.read_file as rf
import numpy as np
import fromage.utils.array_operations as ao
from fromage.utils.mol import Mol
from fromage.utils.atom import Atom

test_dir = os.path.dirname(os.path.abspath(__file__))
test_data_dir = os.path.join(test_dir,'data')

def _in_data(file_name):
    """Return absolute name of file in data/ dir"""
    return os.path.join(test_data_dir,file_name)

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

@pytest.fixture
def newat():
    """Atom object C at origin"""
    return Atom("C", 0.0, 0.0, 0.0)

# Mol fixtures
@pytest.fixture
def at_list():
    """Water dimer atom list"""
    out_list = rf.read_pos(_in_data("h2o_dimer.xyz"))
    return out_list

@pytest.fixture
def h2o_dimer_mol(at_list):
    """Water dimer Mol object"""
    out_mo = Mol(at_list)
    return out_mo

@pytest.fixture
def hc1_mol():
    """HC1 monomer"""
    mol = rf.mol_from_file(_in_data("hc1_mol.xyz"))
    return mol


@pytest.fixture
def hc1_quad():
    """HC1 quadrimer"""
    out_mo = Mol(rf.read_pos(_in_data("hc1_quad.xyz")))
    return out_mo

@pytest.fixture
def hc1_cell():
    """HC1 cell"""
    cell = Mol(rf.read_pos(_in_data("hc1_cell.xyz")))
    vectors = np.array([[12.1199998856, 0.0, 0.0],
                        [0.0, 10.2849998474, 0.0],
                        [-5.4720203118, 0.0, 11.2441994632]])
    cell.vectors = vectors
    return cell

@pytest.fixture
def hc1_complete_cell():
    """HC1 completed cell"""
    cell = Mol(rf.read_pos(_in_data("hc1_complete_cell.xyz")))
    vectors = np.array([[12.1199998856, 0.0, 0.0],
                        [0.0, 10.2849998474, 0.0],
                        [-5.4720203118, 0.0, 11.2441994632]])
    cell.vectors = vectors
    return cell
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
    out_mol = rf.mol_from_file(_in_data("h2o_repeated.xyz"))
    return out_mol

@pytest.fixture
def rectangle_mol():
    """Roughly rectangular soup of atoms"""
    mol = rf.mol_from_file(_in_data("rectangle_atom_soup.xyz"))
    return mol

# Array fixtures
@pytest.fixture
def h2o_dim_array():
    """Coordinate array for H2O dimer"""
    mol = rf.mol_from_file(_in_data("h2o_dimer.xyz"))
    arr = mol.coord_array()
    return arr

@pytest.fixture
def h2o_dim_dist_arr(h2o_dim_array):
    """Distance array for the H2O dimer"""
    return(ao.dist_mat(h2o_dim_array))

@pytest.fixture
def hc1_array():
    """Coordinate array for the HC1 monomer"""
    mol = rf.mol_from_file(_in_data("hc1_mol.xyz"))
    arr = mol.coord_array()
    return arr

@pytest.fixture
def rectangle_array():
    """Coordinate array for corners of a rectangle"""
    lis_coord = [[0.,0.,0.],[0.,2.,0.],[4.,2.,0.],[4.,0.,0.]]
    arr = np.array(lis_coord)
    return arr

@pytest.fixture
def z_plane_coeffs():
    plane_coeffs = np.array([0.,0.,1.,0.])
    return plane_coeffs

@pytest.fixture
def arbitrary_pair():
    """Pair of coordinates"""
    lis_pair = [[-0.5,0.5,3.],[3.,0.5,-1.]]
    arr = np.array(lis_pair)
    return arr

@pytest.fixture
def arbitrary_flat_points():
    """Sea of points in z=0"""
    lis_coord = [[1.,2.,0.],
                [0.1,-0.1,0.],
                [3.,1.,0.],
                [6.1,3.,0.],
                [-0.1,3.1,0.],
                [6.1,-0.1,0.]]
    arr = np.array(lis_coord)
    return arr

@pytest.fixture
def arbitrary_flat_vertices():
    """Vertices of a quadrangle in z=0"""
    lis_vert = [[-0.1,3.1,0.],
                [6.1,3.,0.],
                [6.1,-0.1,0.],
                [0.1,-0.1,0.]]
    arr = np.array(lis_vert)
    return arr

@pytest.fixture
def arbitrary_vertices():
    """Vertices of a quadrangle"""
    lis_vert = [[-0.5,0.5,3.],[0.5,2.,0.],[3.,0.5,-1.],[0.5,-0.5,5.]]
    arr = np.array(lis_vert)
    return arr

@pytest.fixture
def triangle_shape_coord():
    """Set of coordinates in a triangle shape"""
    lis_coord = [[2.,2.,0.],
                 [1.,1.,0.],
                 [1.,4.,0.],
                 [6.,2.,0.],
                 [3.,1.,0.]]
    arr = np.array(lis_coord)
    return arr

@pytest.fixture
def triangle_corners_4():
    lis_coord = [[1.,4.,0.],
                 [6.,2.,0.],
                 [6.,2.,0.],
                 [1.,1.,0.]]
    arr = np.array(lis_coord)
    return arr

@pytest.fixture
def vec_100():
    """The (1,0,0) vector"""
    return np.array([1.,0.,0.])

@pytest.fixture
def vec_220():
    """The (2,2,0) vector"""
    return np.array([2.,2.,0.])

@pytest.fixture
def vec_111():
    """The (1,1,1) vector"""
    return np.array([1.,1.,1.])

# Dimer fixtures
@pytest.fixture
def rectangle_dimer():
    """Dimer object of two rectangle-like molecules"""
    dim = rf.dim_from_file(_in_data("rectangle_mol_dimer.xyz"))
    return dim

@pytest.fixture
def h2o_dimer():
    """Dimer object of two water"""
    dim = rf.dim_from_file(_in_data("h2o_dimer.xyz"))
    return dim
