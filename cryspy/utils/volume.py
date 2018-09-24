import numpy as np

import cryspy.io.edit_file as ef
import cryspy.utils.per_table as pt

class CubeGrid(object):
    """
    A grid of voxels with attached values for each one

    The grid starts at the origin and propagates in a parallelepiped

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

    def get_enclosing_vectors(self):
        n_vox = np.array([self.x_num, self.y_num, self.z_num])
        enclosing_vectors = (self.vectors.T*n_vox).T
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
        """Generate a grid from its centre , box dimension, and resolution

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
                r = atom_i.dist2(*point[:3])
                if scaled:
                    r *= atom_i.vdw**2
                if r < min_dist2:
                    min_dist2 = r
                    close_to_mol = True
            for atom_j in rest:
                r = atom_j.dist2(*point[:3])
                if scaled:
                    r *= atom_i.vdw**2
                if r < min_dist2:
                    min_dist2 = r
                    close_to_mol = False
            if close_to_mol:
                point[3] = 1
            else:
                point[3] = 0

        return

    def vdw_vol(self, mol):
        """Give each point in the grid a value of 1 if it is inside the vdw
        radius of one of the atoms in the molecule"""

        for point in self.grid:
            # empty grid first
            point[3] = 0
            for atom in mol:
                if atom.dist2(*point[:3]) < atom.vdw**2:
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

    def shell_region(self, sample_atoms, inner_r, outer_r):
        """
        Return grid points in shell regions around given points

        The shell region is determined by inner and outer radii which are then
        scaled by the wdv radii of the corresponding atoms.

        Parameters
        ----------
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
        for point in self.grid:
            add = False
            for atom in sample_atoms:
                # we compare squared distances to limit the amount of sqrt operations
                in_r_scaled2 = (inner_r * atom.vdw)**2
                out_r_scaled2 = (outer_r * atom.vdw)**2
                dist2 = atom.dist2(point[0],point[1],point[2])
                if in_r_scaled2 <= dist2 <= out_r_scaled2:
                    add = True
                    break
            if add:
                shell_points.append(point.tolist())
        np.array(shell_points)
        return shell_points

    def supergrid(self):
        """
        Expand the grid so that the new grid has 8 times the volume

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
        null_vec = np.array([0,0,0])
        grids = []
        for trans_a in [null_vec, lattice_vectors[0]]:
            for trans_c in [null_vec, lattice_vectors[1]]:
                for trans_b in [null_vec, lattice_vectors[2]]:
                    new_grid = self.grid.copy()
                    new_grid[:,0:3] -= (trans_a + trans_b + trans_c)
                    grids.append(new_grid)

        unsorted = np.concatenate(grids)
        # The cube values are not yet order like:
        #for i_x in x(for i_y in y(for i_z in z))
        self.grid = unsorted[np.lexsort(np.rot90(unsorted))]
        self.origin = -lattice_vectors.sum(axis=0)
        self.x_num *= 2
        self.y_num *= 2
        self.z_num *= 2

        return self

    def confine_grid(self):
        """Confine grid points in the enclosing box"""
        for row in self.grid:
            pos = row[0:3]
        pass
