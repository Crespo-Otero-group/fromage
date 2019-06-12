import numpy as np

import fromage.io.edit_file as ef
from copy import deepcopy


class CubeGrid(object):
    """
    A grid of voxels with attached values for each one

    This objects contains the information present in a cube file minus the
    atoms inside of it. Its main component is the self.grid but it contains
    additional information which would be compuationally costly to extract
    from the grid: the origin, the lattice vectors, the grid spacing. This
    oject is delicate to use because in order to be writeable as a cube file,
    the self.grid must a) be ordered and b) be on points of the grid. Therefore
    not any translation of the grid will do. Here are provided functions which
    break these rules along with ones which repair them. self.confine_sort is
    therefore important to use on objects which will later be printed

    Attributes
    ----------
    vectors : 3x3 numpy array
        The three vectors defining the shape of one voxel in Angstrom
    x_num, y_num, z_num : ints
        Length of the parallelepiped in units of voxels
    origin : numpy array of length 3
        Origin of the parallelepiped in Angstrom
    grid : numpy array of x_num * y_num * z_num * x 4 dimension
        This determines a value and position at each point in space which is
        the origin of a voxel. The format is [[x1,y1,z1,val],[x1,y1,z2,val],...]

    """

    def __init__(self, vectors=np.zeros((3, 3)), x_num=1, y_num=1, z_num=1, origin=np.array([0.0, 0.0, 0.0])):
        try:
            self.vectors = np.array(vectors)
        except ValueError:
            print("Voxel vectors could not be cast to 3x3 numpy array")

        try:
            self.x_num = int(x_num)
            self.y_num = int(y_num)
            self.z_num = int(z_num)
        except ValueError:
            print("There must be an integer number of voxels")

        try:
            self.origin = np.array(origin)
        except ValueError:
            print("The origin coordinates could not be cast to a numpy array")

        self.dimension = self.x_num * self.y_num * self.z_num
        # initiate grid
        self.grid = np.zeros((self.dimension, 4))

    def copy(self):
        return deepcopy(self)

    def get_enclosing_vectors(self):
        n_vox = np.array([self.x_num, self.y_num, self.z_num])
        enclosing_vectors = (self.vectors.T * n_vox).T
        return enclosing_vectors

    def set_grid_coord(self):
        """Generate the coordinates of the voxel origins on the grid"""

        count = 0
        for x_i in np.arange(0, self.x_num):
            for y_i in np.arange(0, self.y_num):
                for z_i in np.arange(0, self.z_num):
                    entry = np.append(
                        np.dot(self.vectors.T, np.array([x_i, y_i, z_i])) + self.origin, 0)
                    self.grid[count] = entry
                    count += 1
        return

    def grid_from_point(self, x, y, z, res=10, box=np.array([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]])):
        """
        Generate a grid from its centre, box dimension, and resolution

        Parameters
        ----------
        x, y, z : floats
            Centre of the grid
        res : int
            number of voxels per side of the box
        box : 3x3 numpy array
            The three vectors defining the bounding box for the parallelepiped

        """
        # Find centre of box
        cen = (box[0] + box[1] + box[2]) / 2
        self.origin = np.array([x, y, z]) - cen
        self.vectors = box / res
        self.x_num = self.y_num = self.z_num = res
        self.dimension = self.x_num * self.y_num * self.z_num
        self.grid = np.zeros((self.dimension, 4))

        return

    def proximity(self, mol, rest, scaled=True):
        """
        Give each point in the grid a value of 1 if it is closest to the molecule

        As a side effect, the origin is translated by half a voxel so that the
        centre of each voxel corresponds to its value as opposed to the origin of
        each voxel.

        Parameters
        ----------
        mol: list of Atom objects
            The central molecule which we want to enclose in the grid
        rest: list of Atom objects
            The rest of the atoms

        """
        for point in self.grid:
            close_to_mol = False
            min_dist2 = float("inf")
            for atom_i in mol:
                r = atom_i.c_dist2(*point[:3])
                if scaled:
                    r /= atom_i.vdw**2
                if r < min_dist2:
                    min_dist2 = r
                    close_to_mol = True
            for atom_j in rest:
                r = atom_j.c_dist2(*point[:3])
                if scaled:
                    r /= atom_j.vdw**2
                if r < min_dist2:
                    min_dist2 = r
                    close_to_mol = False
                    break
            if close_to_mol:
                point[3] = 1
            else:
                point[3] = 0

        return

    def vdw_vol(self, mol):
        """Give each point in the grid a value of 1 if it is inside the vdw radius of one of the atoms in the molecule"""

        for point in self.grid:
            # empty grid first
            point[3] = 0
            for atom in mol:
                if atom.c_dist2(*point[:3]) < atom.vdw**2:
                    point[3] = 1
                    break
        return

    def subtract_grid(self, in_grid):
        """
        Remove the results of another grid from the current grid.

        Parameters
        ----------
        in_grid : numpy N x 4 array
            Same dimensions as self.grid

        """

        for i, j in zip(self.grid, in_grid):
            i[3] -= j[3]
        return

    def add_grid(self, in_grid):
        """
        Add the results of another grid to the current grid.

        Parameters
        ----------
        in_grid : numpy N x 4 array
            Same dimensions as self.grid

        """

        for i, j in zip(self.grid, in_grid):
            i[3] += j[3]
        return

    def out_cube(self, file_name, atoms):
        """Write a cube file with the current state of the grid"""
        values = np.array([point[3] for point in self.grid])
        ef.write_cube(file_name, self.origin, self.vectors, self.x_num,
                      self.y_num, self.z_num, atoms, values)
        return

    def volume(self):
        filled = 0
        for entry in self.grid:
            if entry[3] != 0:
                filled += 1

        vox_vol = np.linalg.det(self.vectors)

        return filled * vox_vol

    def expand(self):
        """
        BROKEN! Expand the grid so that the new grid has 8 times the volume

        The original grid has the origin at 0,0,0 so we add grids at origins:
        -a, 0, 0
        0,-b, 0
        0, 0,-c
        -a,-b, 0
        -a, 0,-c
        0,-b,-c
        -a,-b,-c

        """
        lattice_vectors = self.get_enclosing_vectors()
        null_vec = np.array([0, 0, 0])
        grids = []
        for trans_a in [null_vec, lattice_vectors[0]]:
            for trans_c in [null_vec, lattice_vectors[1]]:
                for trans_b in [null_vec, lattice_vectors[2]]:
                    new_grid = self.grid.copy()
                    new_grid[:, 0:3] -= (trans_a + trans_b + trans_c)
                    grids.append(new_grid)

        unsorted = np.concatenate(grids)
        # The cube values are not yet order like:
        # for i_x in x(for i_y in y(for i_z in z))
        self.grid = unsorted[np.lexsort(np.rot90(unsorted))]
        self.origin = -lattice_vectors.sum(axis=0)
        self.x_num *= 2
        self.y_num *= 2
        self.z_num *= 2

        return self

    def supergrid_unsorted(self, trans):
        """
        Make a supercell cube out of the original cube grid but no sorting

        No sorting means that the values will not be in the correct order for
        the cube file.

        Parameters
        ----------
        trans : numpy array of length 3
            Multiplications of the primitive cell

        """
        grids = []
        lattice_vectors = self.get_enclosing_vectors()
        for a_mult in range(trans[0]):
            for b_mult in range(trans[1]):
                for c_mult in range(trans[2]):
                    new_grid = self.grid.copy()
                    new_grid[:, 0:3] += a_mult * lattice_vectors[0] + \
                        b_mult * lattice_vectors[1] + \
                        c_mult * lattice_vectors[2]
                    grids.append(new_grid)

        # The cube values are not yet order like:
        # for i_x in x(for i_y in y(for i_z in z))
        self.grid = np.concatenate(grids)
        self.x_num *= trans[0]
        self.y_num *= trans[1]
        self.z_num *= trans[2]

        return self

    def supergrid(self, trans):
        """
        Make a supercell cube out of the original cube grid with sorting

        This is the function to call if the grid is then to be made into a cube
        file.

        Parameters
        ----------
        trans : numpy array of length 3
            Multiplications of the primitive cell

        """
        self.supergrid_unsorted(trans)
        self.confine_sort()

        return self

    def dir_to_frac_pos(self):
        """
        Move all grid points to fractional coordinates

        This breaks the assumption that the grid is in real space. Therefore
        this function is to be used with care
        """
        new_grid = self.grid.copy()
        #new_grid[:, 0:3] -= self.origin
        lattice_vectors = self.get_enclosing_vectors()
        # transpose to get the transformation matrix
        M = np.transpose(lattice_vectors)
        # inverse transformation matrix
        U = np.linalg.inv(M)
        # get a matrix A such that A[i] = U dot self.grid[i] excluding the 4th
        # column which remains intact.
        new_grid[:, 0:3] = np.einsum('ij,kj->ki', U, self.grid[:, 0:3])
        self.grid = new_grid
        return

    def sort_adjust_frac_pos(self, sorting=True, rounding=True):
        """
        Put the fractional grid on fractional grid points and sort

        The idea here is that in order to be able to use cube files, the grid
        needs to be sorted in a specific way and the real space points must
        coincide with grid points of the cube file, otherwise hard to diagnose
        bugs happen. However in order to match grid points, slight translations
        need to be done in the form of rounding, thus making the grid points
        inaccurate for use in computation. Therefore this step is only to be
        used when the end result is a cube file, not a grid for computation.

        """
        new_grid = self.grid.copy()
        # now we make sure that the points are on grid points of the mesh
        xyz_nums = np.array([self.x_num, self.y_num, self.z_num, ])
        new_grid[:, 0:3] *= xyz_nums
        if rounding:
            new_grid[:, 0:3] = np.round(new_grid[:, 0:3], 0)
        new_grid[:, 0:3] /= xyz_nums
        # And now confine them to the cell
        new_grid[:, 0:3] = np.mod(new_grid[:, 0:3], 1)
        # get in the proper order for cube files
        if sorting:
            new_grid = new_grid[np.lexsort(np.rot90(new_grid))]
        self.grid = new_grid
        return

    def frac_to_dir_pos(self):
        """Move all grid points to direct coordinates"""
        lattice_vectors = self.get_enclosing_vectors()  # + self.origin
        # print(self.origin)
        # see dir_to_frac_pos for detils
        new_grid = np.einsum('ij,ki->kj', lattice_vectors, self.grid[:, 0:3])
        #self.grid[:, 0:3] = new_grid + self.origin
        self.grid[:, 0:3] = new_grid
        return

    def confine_sort(self):
        """
        Confine grid points in the enclosing box in the correct order

        The enclosing box has origin in self.origin and bounding vectors
        self.get_enclosing_vectors().

        """
        self.dir_to_frac_pos()
        self.sort_adjust_frac_pos()
        self.frac_to_dir_pos()
        return

    def confine_unordered(self):
        """
        Confine grid points in the enclosing box but lose the correct order

        This means that the grid is no longer ordered or on grid points defined
        by the voxel spacing of the CubeGrid.
        """
        self.dir_to_frac_pos()
        self.sort_adjust_frac_pos(sorting=False, rounding=False)
        self.frac_to_dir_pos()

    def translate_grid(self, trans_vec):
        """Translate grid and origin by a numpy vector"""
        self.grid[:, 0:3] += trans_vec
        #self.origin += trans_vec
        return

    def translate_inplace(self, trans_vec):
        """
        Translate the cell and then confine it back to the enclosing box

        This allows for the grid to move in the periodic cell without having to
        change the origin

        Parameters
        ----------
        trans_vec : 3x1 numpy array
            The vector by which to translate the grid

        """
        self.translate_grid(trans_vec)
        self.confine_sort()
        # self.translate_grid(-trans_vec)
        return

    def unord_trans_inplace_grid(self, trans_vec):
        """Get the unordered Cub after an inplace translation"""
        fresh_cub = self.copy()
        fresh_cub.translate_grid(trans_vec)
        fresh_cub.confine_unordered()
        return fresh_cub

    def centered_quad(self, trans_vec):
        """
        Produce a 4x4x4 supercell centered at the origin after inplace translate

        The point is to have a large grid which encloses a specific point of the
        original grid. This way we can center a molecule (or several) at the
        origin and sample as much as we want without risking hitting the
        supercell walls

        Parameters
        ----------
        trans_vec : 3x1 numpy array
            The vector by which to translate inplace. The point which was
            originally at trans_vec ends up at the origin

        """
        new_cub = self.copy()
        new_cub.translate_inplace(trans_vec)
        super_cub = new_cub.supergrid([4, 4, 4])

        center = np.sum(super_cub.get_enclosing_vectors(), axis=0) / 2
        super_cub.origin -= center

        return super_cub
