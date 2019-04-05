import numpy as np
import fromage.utils.array_operations as ao

class GeomInfo(object):
    """
    Contains the extended geometrical information for Mol objects

    This object is meant exclusively to store geometrical data. It needs to be
    initialised as empty since calculating some of the properties is time
    intensive. As such, it needs to be assignable which rules out simply using
    a namedtuple.

    Attributes
    ----------
    coord_array : np array of shape Nat X 3
        Atom coordinates
    plane_coeffs : 3 x 1 np array
        Coefficients for the plane which averages coord_array
    prin_ax : 3 x 1 np array
        Vector representing the principal axis of the molecule
    sec_ax : 3 x 1 np array
        Vector representing the secondary axis of the molecule
    perp_ax : 3 x 1 np array
        Vector perpendicular to the other two such that
        perp_ax = prin_ax (cross) sec_ax
    ignore_kinds : list of Atom.kinds
        Specific atoms to be ignored in the geometry parameters
    ignore_hydrogens : bool
        Ignore hydrogens when calculating geometry parameters

    """
    def __init__(self):
        self.coord_array = np.array([0])
        self.plane_coeffs = np.array([0])
        self.prin_ax = np.array([0])
        self.sec_ax = np.array([0])
        self.perp_ax = np.array([0])
        self.ignore_kinds = []
        self.ignore_hydrogens = False
        # this is useful for determining axes differently
        self.linear = False

    def __str__(self):
        out_str = "Coordinate array:\n" + str(self.coord_array) + "\nPlane coefficients:\n" + str(
            self.plane_coeffs) + "\nPrincipal axis:\n" + str(self.prin_ax) + "\nSecondary axis:\n" + str(
            self.sec_ax) + "\nPerpendicular axis:\n" + str(self.perp_ax)
        return out_str

    def __repr__(self):
        return self.__str__()


def coord_array(self):
    """
    Return a numpy array of the coordinates

    Returns
    -------
    coord_arr : Nat x 3 numpy array
        Array of the form [[x1,y1,z1],[x2,y2,z2],...]

    """
    if self.geom.ignore_kinds:
        self.set_connectivity
    list_coord = []
    nat = len(self)
    coord_arr = np.zeros((nat,3))
    for atom in self:
        add_atom = True
        # potentially remove hydrogen
        if self.geom.ignore_hydrogens:
            if atom.elem == 'H':
                add_atom = False
        # potentially remove a kind of atom
        if self.geom.ignore_kinds:
            if atom.kind in self.geom.ignore_kinds:
                add_atom = False
        if add_atom:
            new_row = [atom.x,atom.y,atom.z]
            list_coord.append(new_row)

    coord_arr = np.array(list_coord)
    return coord_arr

def calc_coord_array(self):
    """Set the coordinate array in geom"""
    self.geom.coord_array = self.coord_array()

def plane_coeffs(self):
    """
    Return numpy array of the plane coefficients which average the coords

    Returns
    -------
    plane_coeffs : length 4 numpy array
        Plane equation coefficients such that a point on the plane is:
        ax + by + cz + d = 0. The array is [a,b,c,d]

    """
    if np.count_nonzero(self.geom.coord_array) == 0:
        self.calc_coord_array()
    plane_coeffs = ao.plane_from_coord(self.geom.coord_array)
    return plane_coeffs

def calc_plane_coeffs(self):
    """Set the plane coefficients in geom"""
    self.geom.plane_coeffs = self.plane_coeffs()

def axes(self):
    """
    Return principal, secondary and perpendicular axes of the Mol

    Returns
    -------
    axes_out : 3 x 3 np array
        Principal, followed by secondary and perpendicular axes of the molecule

    """
    if np.count_nonzero(self.geom.plane_coeffs) == 0:
        self.calc_plane_coeffs()
    # get the quadrangle which best describes the coordinates (possibly a
    # triangle with the far point repeated twice)
    vertices = ao.quadrangle_from_coord(self.geom.coord_array)
    # if the mol is linear, we need to reorder the quadrangle
    if self.geom.linear:
        vertices = np.array([vertices[1],vertices[2],vertices[3],vertices[0]])
    # if the mole is rectangular, we want an embedded quadrangle
    else:
        # get the embedded quadrangle vertices
        vertices = ao.embedded_vert(vertices)
    # get vectors from projected diagonals
    axes_out_raw = ao.project_quad_to_vectors(vertices,self.geom.plane_coeffs)
    # we want the first raw to be the secondary and vice versa and the principal
    # to be *(-1) in order to maintain a convention
    axes_out_unnormal = np.array([-axes_out_raw[1], axes_out_raw[0]])
    # orthonogalise them
    # if the mol is linear, just move the secondary axis
    if self.geom.linear:
        axes_out_prin_sec = ao.orthogonalise_asym(axes_out_unnormal)
    # if the mol is not linear, move principal and secondary axes equally
    else:
        axes_out_prin_sec = ao.orthogonalise_sym(axes_out_unnormal)
    # get the perpendicular vector
    perp = np.cross(axes_out_prin_sec[0], axes_out_prin_sec[1])
    # ensure normalisation
    perp = perp/np.linalg.norm(perp)

    lis_axes_out = list(axes_out_prin_sec)
    lis_axes_out.append(perp)
    axes_out = np.array(lis_axes_out)
    return axes_out

def calc_axes(self):
    """Set the principal and secondary axes in geom"""
    axes = self.axes()
    self.geom.prin_ax = axes[0]
    self.geom.sec_ax = axes[1]
    self.geom.perp_ax = axes[2]
