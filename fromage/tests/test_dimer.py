import fromage.io.read_file as rf
import numpy as np
from fromage.utils.dimer import Dimer
from fromage.tests.conftest import _in_data
from numpy.testing import assert_allclose
from pytest import approx


def test_init(h2o_dimer_mol):
    mol_a = h2o_dimer_mol.select(0)
    mol_b = h2o_dimer_mol.select(2)
    dim = Dimer(mol_a, mol_b)
    assert len(dim.mol_a) == 3
    assert len(dim.mol_b) == 3


def test_dimer_from_file():
    dim = rf.dimer_from_file(_in_data("h2o_dimer.xyz"))
    assert len(dim.mol_a) == 3
    assert len(dim.mol_b) == 3


def test_rectangle_dim_angles(rectangle_dimer):
    rectangle_dimer.calc_angles()
    assert rectangle_dimer.alpha == approx(18.60637, rel=10e-4)
    assert rectangle_dimer.beta == approx(29.3625, rel=10e-4)
    assert rectangle_dimer.gamma == approx(34.2641, rel=10e-4)


def test_rectangle_dim_slip(rectangle_dimer):
    ang = rectangle_dimer.slip_angle()
    assert ang == approx(18.036836, rel=10e-4)


def test_h2_dimer_inter_distance(h2_dimer):
    assert h2_dimer.inter_distance() == approx(5.0)
    assert h2_dimer.inter_distance(mode="cov") == approx(4.54)
    assert h2_dimer.inter_distance(mode="vdw") == approx(2.82)
    assert h2_dimer.inter_distance(method="centroid") == approx(6.0)


def test_h2o_dimer_inter_distance(h2o_dimer):
    assert h2o_dimer.inter_distance() == approx(2.3512055, rel=10e-4)
    assert h2o_dimer.inter_distance(mode="cov") == approx(1.4412055, rel=10e-4)
    assert h2o_dimer.inter_distance(mode="vdw") == approx(-0.2587945, rel=10e-4)
    assert h2o_dimer.inter_distance(method="centroid") == approx(3.0413813, rel=10e-4)


def test_he_images(he_dimer):
    vectors = np.array([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]])
    lis = he_dimer.images(vectors)
    assert len(lis) == 27


def test_identical_to(h2o_dimer, h2o_dimer_jumbled):
    assert h2o_dimer.identical_to(h2o_dimer_jumbled)


def test_inter_distance(h2_dimer):
    res = h2_dimer.sorted_inter_distances()
    expected = np.array([5.0, 6.0, 6.0, 7.0])
    assert_allclose(res, expected)


def test_same_geom(h2o_dimer, h2o_dimer_jumbled_trans):
    assert h2o_dimer.same_geom(h2o_dimer_jumbled_trans)
