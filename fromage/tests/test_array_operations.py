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


def test_dist_vectors_three_by_three(three_points):
    """Test the distance vectors between three points and themselves"""
    test_vecs = np.array(
        [
            [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 4.0, 0.0]],
            [[0.0, 0.0, -3.0], [0.0, 0.0, 0.0], [0.0, 4.0, -3.0]],
            [[0.0, -4.0, 0.0], [0.0, -4.0, 3.0], [0.0, 0.0, 0.0]],
        ]
    )

    assert ao.dist_vec(three_points, three_points) == approx(test_vecs)


def test_dist_vectors_three_by_two(three_points, two_points):
    """Test the distance vectors between three and two points"""
    test_vecs = np.array(
        [
            [[0.0, 1.0, 0.0], [0.0, 0.0, 2.0]],
            [[0.0, 1.0, -3.0], [0.0, 0.0, -1.0]],
            [[0.0, -3.0, 0.0], [0.0, -4.0, 2.0]],
        ]
    )

    assert ao.dist_vec(three_points, two_points) == approx(test_vecs)


def test_dist_three_by_three(three_points):
    """Test the distance vectors between three and themselves"""
    test_dists = np.array([[0.0, 3.0, 4.0], [3.0, 0.0, 5.0], [4.0, 5.0, 0.0]])

    assert ao.per_dist_mat(three_points, three_points) == approx(test_dists)


def test_dist_three_by_three_per(three_points, lattice_vectors_555):
    """Test the distance vectors between three and themselves in a cell"""
    test_dists = np.array(
        [[0.0, 2.0, 1.0], [2.0, 0.0, np.sqrt(5)], [1.0, np.sqrt(5), 0.0]]
    )

    assert ao.per_dist_mat(
        three_points, three_points, lat_vec=lattice_vectors_555
    ) == approx(test_dists)

def test_dir_to_frac(vec_100, lattice_vectors_444):
    test_frac = np.array([0.25, 0.0, 0.0])
    frac_coord = ao.dir_to_frac(vec_100, lattice_vectors_444)
    assert_allclose(frac_coord, test_frac)

def test_frac_to_dir(vec_100, lattice_vectors_444):
    frac_coord = ao.dir_to_frac(vec_100, lattice_vectors_444)
    dir_coord = ao.frac_to_dir(frac_coord, lattice_vectors_444)
    assert_allclose(dir_coord, vec_100)

def test_confine(lattice_vectors_444):
    out_vec = np.array([6.0, 0.0, 0.0])
    in_vec = ao.confine(out_vec, lattice_vectors_444)
    test_in_vec = np.array([2.0, 0.0, 0.0])
    assert_allclose(in_vec, test_in_vec)

def test_per_translate(two_points, lattice_vectors_555):
    translation = np.array([0.0, 6.0, 0.0])
    new_points = ao.per_translate(two_points, translation, lattice_vectors_555)
    test_translated = np.array([[0.0, 2.0, 0.0],[0.0, 1.0, 2.0]])
    assert_allclose(new_points, test_translated)


