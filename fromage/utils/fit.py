"""Fit point charges to match a given potential"""
import numpy as np
from fromage.utils.mol import Mol

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
        The sampling points with rows as x1 y1 z1 value1

    """
    shell_points = []
    for point in in_grid:
        add = False
        for atom in sample_atoms:
            # we compare squared distances to limit the amount of sqrt
            # operations
            in_r_scaled2 = (inner_r * atom.cov)**2
            out_r_scaled2 = (outer_r * atom.cov)**2
            dist2 = atom.v_dist2(point[0:3])
            if in_r_scaled2 <= dist2 <= out_r_scaled2:
                add = True
                break
        if add:
            shell_points.append(point.tolist())
    shell_points = np.array(shell_points)
    return shell_points

def alt_shell_region(in_grid, sample_atoms, inner_r, outer_r):
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
        The sampling points with rows as x1 y1 z1 value1

    """
    # array full of False
    keep_bool_arr = np.zeros((1,len(in_grid)),dtype=bool)[0]
    for atom in sample_atoms:
        scaled_inner_r2 = (atom.vdw * inner_r)**2
        scaled_outer_r2 = (atom.vdw * outer_r)**2
        # array of vectors atom -> points
        diff = in_grid[:,0:3] - atom.get_pos()
        # array of distances**2 atom -> points
        distances2 = np.einsum('ij,ij->i',diff,diff)
        bool_arr = (scaled_inner_r2 < distances2) & (distances2 < scaled_outer_r2)
        # make True all of the indices that satisfied the condition
        keep_bool_arr += bool_arr

    shell_points = in_grid[np.where(keep_bool_arr)[0]]
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
#        entry = 1 / point.dist(sample)
        entry = 1 / point.v_dist(sample)

        l_row.append(entry)
    row = np.array(l_row)
    return row


def dep_var(var_points, fix_points, samples):
    """Return the dependent variable array"""
    l_dep = []
    for sam in samples:
        # entry = sam.es - \
        #     var_points.es_pot(sam.get_pos()) - fix_points.es_pot(sam.get_pos())
        entry = sam[3] - \
            var_points.es_pot(sam[0:3]) - fix_points.es_pot(sam[0:3])
        l_dep.append(entry)
    out_dep = np.array(l_dep)
    return out_dep


def fit_points(var_points, samples, fix_points= None):
    """
    Return a new set of point charges that matches the potential at points

    Parameters
    ----------
    var_points : Mol object
        Point charges to be fitted
    samples : list of Mol objects or just one Mol object
        Sampling points each with an associated electrostatic potential
    fix_points : Mol object
        Point charges to remain in place
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

    print("RMSD: "+str(np.sqrt(np.sum(res[3])/len(samples))))
    print(res[1:])
    fitting = res[0]
    print(fitting[0:6])
    var_points.change_charges(var_points.charges() + fitting)

    return var_points

def shells_from_cell(cell_cub, central_mol, trans_vec, inner_r, outer_r):
    """
    Return the points in shell regions of a cube file after translation

    Parameters
    ----------
    cell_cub : CubeGrid object
        The cube file of the unit cell before any translations
    central_mol : Mol object
        The atoms which need to be surrounded by shells. These have already been
        selected from lattice positions and translated by trans_vec
    trans_vec : 3x1 numpy object
        Vector by which cental_mol was translated and which therefore needs to
        translate the cell_cub.grid
    inner_r : float
        The inner radius of the shell region for sampling around atoms of the
        central_mol. This distance is then scaled by wdv radius
    outer_r : float
        The outer radius of the shell region for sampling around atoms of the
        central_mol. This distance is then scaled by wdv radius
    Returns
    -------
    shell_points : numpy Nx4 array
        The sampling points with rows as x1 y1 z1 value1

    """
    trans_cub = cell_cub.unord_trans_inplace_grid(trans_vec)
    sup_cube = trans_cub.supergrid_unsorted([4,4,4])
    center = np.sum(sup_cube.get_enclosing_vectors(),axis=0)/2
    grid = sup_cube.grid
    grid[:,0:3] -= center
    sample_points = alt_shell_region(grid, central_mol, inner_r, outer_r)
    return sample_points

#def fit_clust(in_cell, in_label, inner_r, outer_r):
def fit_clust(in_cell, in_labels, in_cube):
    print("Start test")
    print("Complete mols: start")
    mol, mod_cell, trans = in_cell.centered_mols(in_labels, return_trans = True)
    print("Complete mols: done")

    print("Make cluster: start")
    shell = mod_cell.make_cluster(15)
    print("Make cluster: done")

    print("Remove kernel: start")
    for atom in mol:
        if atom in shell:
            shell.remove(atom)

    print("Remove kernel: done")

    print("Shell sampling: start")
    samples = shells_from_cell(in_cube, mol, trans, 0.9995, 1.0)
    print("Shell sampling: done")

    print("Fitting: start")
    fit_points(shell, samples)
    print("Fitting: done")


    return
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
