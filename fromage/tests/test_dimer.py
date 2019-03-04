import fromage.io.read_file as rf
from fromage.utils.dimer import Dimer
from fromage.tests.conftest import _in_data
from pytest import approx

def test_init(h2o_dimer_mol):
    mol_a = h2o_dimer_mol.select(0)
    mol_b = h2o_dimer_mol.select(2)
    dim = Dimer([mol_a, mol_b])
    assert len(dim.mols[0]) == 3
    assert len(dim.mols[1]) == 3

def test_dim_from_file():
    dim = rf.dim_from_file(_in_data("h2o_dimer.xyz"))
    assert len(dim.mols[0]) == 3
    assert len(dim.mols[1]) == 3

def test_rectangle_dim_angles(rectangle_dimer):
    rectangle_dimer.calc_angles()
    assert rectangle_dimer.alpha == approx(5.187172,rel=10e-4)
    assert rectangle_dimer.beta == approx(5.643125,rel=10e-4)
    assert rectangle_dimer.gamma == approx(3.247717,rel=10e-4)
