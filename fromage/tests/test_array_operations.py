import fromage.utils.array_operations as ao
import numpy as np
from numpy.testing import assert_allclose
from pytest import approx


def test_distance(vec_100, vec_220):
    distance = ao.distance(vec_100, vec_220)
    assert distance == approx(2.236, rel=0.001)


def test_closest(vec_100, vec_220):
    origin = np.array([0.0, 0.0, 0.0])
    close_point, dis = ao.closest(origin, [vec_100, vec_220])
    assert_allclose(close_point, vec_100)
    assert dis == approx(1.0)


def test_dist_mat(h2o_dim_dist_arr):
    dis = ao.dist_mat(h2o_dim_dist_arr)
    assert dis[1][0] == approx(0.91, rel=0.1)


def test_find_largest(h2o_dim_dist_arr):
    largest_inds = ao.find_largest(h2o_dim_dist_arr, 3)
    assert largest_inds == [(4, 0), (4, 2), (5, 0)]


# Angle tests
def test_angle_between_vectors(vec_100, vec_220):
    ang = ao.vec_angle(vec_100, vec_220)
    assert ang == approx(45.0)


def test_angle_between_vectors_rad(vec_100, vec_220):
    ang = ao.vec_angle(vec_100, vec_220, degrees=False)
    assert ang == approx(45.0 * np.pi / 180)


# Orthog test
def test_orthogonalise_sym():
    vectors = np.array([[0.1, 0.0, 1.0], [-0.1, 0.0, 1.0]])

    o_vecs = ao.orthogonalise_sym(vectors)
    ovec_1 = o_vecs[0] / np.linalg.norm(o_vecs[0])
    ovec_2 = o_vecs[1] / np.linalg.norm(o_vecs[1])

    expected_1 = np.array([np.sqrt(2) / 2, 0.0, np.sqrt(2) / 2])
    expected_2 = np.array([-np.sqrt(2) / 2, 0.0, np.sqrt(2) / 2])

    assert_allclose(ovec_1, expected_1)
    assert_allclose(ovec_2, expected_2)


def test_orthogonalise_asym():
    vectors = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 1.0]])

    o_vecs = ao.orthogonalise_asym(vectors)
    ovec_1 = o_vecs[0] / np.linalg.norm(o_vecs[0])
    ovec_2 = o_vecs[1] / np.linalg.norm(o_vecs[1])

    expected_1 = np.array([0.0, 0.0, 1.0])
    expected_2 = np.array([0.0, 1.0, 0.0])

    assert_allclose(ovec_1, expected_1)
    assert_allclose(ovec_2, expected_2, atol=1e-05)


def test_rmsd():
    ref = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    alt = np.array([1.0, -1.0, -1.0, -1.0, 1.0])
    assert ao.rmsd(ref, alt) == 1.0


def test_coord_rmsd():
    """Test the rmsd of two points"""
    ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    alt = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    assert ao.coord_rmsd(ref, alt) == approx(np.sqrt(0.5))


def test_dist_vectors():
    """Test the zmatrix between three points and themselves"""
    points = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 4.0, 0.0]])
    test_zmat = np.array([[0.0, 3.0, 4.0], [3.0, 0.0, 5,], [4.0, 5.0, 0.0]])

    ao.dist_vec(points, points) == approx(test_zmat)
