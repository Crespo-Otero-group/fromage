import numpy as np

def cross_product_matrix(in_vector):
    """
    Return the cross product matrix of a 3d vector

    For more information, see:
    https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication

    Parameters
    ----------
    in_vector : 3 x 1 numpy matrix

    """
    I = np.identity(3)
    out_mat = np.zeros((3,3))
    for row in I:
        cross = np.cross(in_vector,row)
        outer = np.outer(cross,row)
        out_mat += outer
    return out_mat


def rotation_matrix(axis_vector, angle, degrees = True):
    """
    Return the rotation matrix corresponding to a rotation axis and angle

    For more information, see:
    https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    Parameters
    ----------
    axis_vector : 3 x 1 numpy array
        A unit vector of the axis of rotation
    angle : float
        Angle of rotation in degrees unless otherwise specified
    degrees : bool (optional)
        Choose between units of degrees of radians. Default True so degrees

    """
    ang = angle
    if degrees:
        ang = np.degrees(ang)
    ux, uy, uz = axis_vector
    return
