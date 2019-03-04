import fromage.io.read_file as rf
import numpy as np
from fromage.utils.dimer import Dimer
from fromage.tests.conftest import _in_data
from pytest import approx

def test_init(h2o_dimer_mol):
    mol_a = h2o_dimer_mol.select(0)
    mol_b = h2o_dimer_mol.select(2)
    dim = Dimer(mol_a, mol_b)
    assert len(dim.mol_a) == 3
    assert len(dim.mol_b) == 3

def test_dim_from_file():
    dim = rf.dim_from_file(_in_data("h2o_dimer.xyz"))
    assert len(dim.mol_a) == 3
    assert len(dim.mol_b) == 3

def test_rectangle_dim_angles(rectangle_dimer):
    rectangle_dimer.calc_angles()
    assert rectangle_dimer.alpha == approx(5.187172,rel=10e-4)
    assert rectangle_dimer.beta == approx(5.643125,rel=10e-4)
    assert rectangle_dimer.gamma == approx(3.247717,rel=10e-4)

def test_h2_dimer_inter_distance(h2_dimer):
    assert h2_dimer.inter_distance() == approx(5.)
    assert h2_dimer.inter_distance(mode='cov') == approx(4.54)
    assert h2_dimer.inter_distance(mode='vdw') == approx(2.82)
    assert h2_dimer.inter_distance(method='centroid') == approx(6.)

def test_h2o_dimer_inter_distance(h2o_dimer):
    assert h2o_dimer.inter_distance() == approx(2.3512055,rel=10e-4)
    assert h2o_dimer.inter_distance(mode='cov') == approx(1.4412055,rel=10e-4)
    assert h2o_dimer.inter_distance(mode='vdw') == approx(-0.2587945,rel=10e-4)
    assert h2o_dimer.inter_distance(method='centroid') == approx(3.0413813,rel=10e-4)

def test_he_images(he_dimer):
    vectors = np.array([[20.0,0.0,0.0],
                        [0.0,20.0,0.0],
                        [0.0,0.0,20.0]])
    lis = he_dimer.images(vectors)
    print(lis)
