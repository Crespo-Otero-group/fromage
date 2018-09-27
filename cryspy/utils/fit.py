"""Fit point charges to match a given potential"""
import numpy as np
from cryspy.utils.mol import Mol

def shell_region(in_grid, sample_atoms, inner_r, outer_r):
    """
    Return grid points in shell regions around given points

    The shell region is determined by inner and outer radii which are then
    scaled by the wdv radii of the corresponding atoms.

    Parameters
    ----------
    in_grid : numpy array of N x 4
        A real space grid representing a field with each row being
        [x y z value]
    sample_atoms : Mol object
        The atoms which are to be enclosed by the shells
    inner_r : float
        The inner radius of the shell before wdv scaling
    outer_r : float
        The outer radius of the shell before scaling
    Returns
    -------
    shell_points : numpy N x 4 array

    """
    shell_points = []
    for point in in_grid:
        add = False
        for atom in sample_atoms:
            # we compare squared distances to limit the amount of sqrt
            # operations
            in_r_scaled2 = (inner_r * atom.vdw)**2
            out_r_scaled2 = (outer_r * atom.vdw)**2
            dist2 = atom.dist2(point[0], point[1], point[2])
            if in_r_scaled2 <= dist2 <= out_r_scaled2:
                add = True
                break
        if add:
            shell_points.append(point.tolist())
    np.array(shell_points)
    return shell_points

def coeff_mat(var_points, samples):
    """Return the coefficients matrix"""
    out_mat = np.zeros((len(samples), len(var_points)))
    for i, sam in enumerate(samples):
        out_mat[i] = coeff_row(var_points, sam)
    return out_mat


def coeff_row(var_points, sample):
    """Return row of the coefficients matrix"""
    l_row = []
    for point in var_points:
        entry = 1 / point.dist_at(sample)
        l_row.append(entry)
    row = np.array(l_row)
    return row


def dep_var(var_points, fix_points, samples):
    """Return the dependent variable array"""
    l_dep = []
    for sam in samples:
        entry = sam.es - \
            var_points.es_pot(sam.get_pos()) - fix_points.es_pot(sam.get_pos())
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
    samples : list of Mol objects or just one Mol object
        Sampling points each with an associated electrostatic potential
    Returns
    -------
    out_points : Mol object
        The points at their same position but with optimised charge value

    """
    if fix_points is None:
        fix_points = Mol([])
    coeffs = coeff_mat(var_points, samples)
    deps = dep_var(var_points, fix_points, samples)

    res = np.linalg.lstsq(coeffs, deps, rcond=None)
    # print(res)
    fitting = res[0]
    var_points.change_charges(var_points.charges() + fitting)

    return var_points


    # atoms = rf.read_pos(cell_file)
    # output_file.write("Read " + str(len(atoms)) + " atoms in cell_file\n")
    # output_file.flush()
    # # the molecule of interest and the atoms which now contain
    # # the full, unchopped molecule
    # # NB: all objects in mol are also referenced inside atoms
    # mol, atoms = ha.complete_mol(max_bl, atoms, atom_label, vectors)
    #
    # # find the centroid of the molecule
    # c_x, c_y, c_z = ha.find_centroid(mol)
    # # translate the molecule and atoms to the centroid
    # for atom in atoms:
    #     atom.translate(-c_x, -c_y, -c_z)
    #
    # # write useful xyz and new cell
    # ef.write_xyz("mol.init.xyz", mol)
    # if print_tweak:
    #     ef.write_xyz("tweaked_cell.xyz", atoms)
