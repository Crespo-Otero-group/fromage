"""Fit point charges to match a given potential"""
import numpy as np

def coeff_mat(points,samples):
    """Return the coefficients matrix"""
    out_mat = np.zeros((len(samples),len(points)))
    for i, sam in enumerate(samples):
        out_mat[i] = coeff_row(points, sam)
    return out_mat

def coeff_row(points, sample):
    """Return row of the coefficients matrix"""
    l_row = []
    for point in points:
        entry = 1/point.dist_at(sample)
        l_row.append(entry)
    row = np.array(l_row)
    return row

def dep_var(points, samples):
    """Return the dependent variable array"""
    l_dep = []
    for sam in samples:
        entry = sam.es - points.es_pot(sam.get_pos())
        l_dep.append(entry)
    out_dep = np.array(l_dep)
    return out_dep

def fit_points(var_points, fix_points, samples):
    """
    Return a new set of point charges that matches the potential at points

    Parameters
    ----------
    points : Mol object
        Points charges to be fitted
    samples : Mol object
        Sampling points each with an associated electrostatic potential
    Returns
    -------
    out_points : Mol object
        The points at their same position but with optimised charge value

    """
    coeffs = coeff_mat(var_points,samples)
    deps = dep_var(var_points,samples)

    res = np.linalg.lstsq(coeffs, deps, rcond=None)[0]
    var_points.change_charges(var_points.charges()+res)
    print(var_points.es_pot([0,0,0]))
    return var_points
