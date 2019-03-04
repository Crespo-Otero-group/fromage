import numpy as np

def cross_product_matrix(in_vector):
    """
    Return the cross product matrix of a 3d vector

    For more information, see:
    https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication

    Parameters
    ----------
    in_vector : 3 x 1 numpy matrix
        Input vector
    Returns
    -------
    out_mat = 3 x 3 numpy array
        The cross product matrix of in_vector

    """
    # Alternate implementation for curiosity
    # I = np.identity(3)
    # out_mat = np.zeros((3,3))
    # for row in I:
    #     cross = np.cross(in_vector,row)
    #     outer = np.outer(cross,row)
    #     out_mat += outer
    out_mat = np.array([[0.,-in_vector[2],in_vector[1]],
                        [in_vector[2],0,-in_vector[0]],
                        [-in_vector[1],in_vector[0],0.]])
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
    Returns
    -------
    rot_mat : 3 x 3 numpy array
        Rotation matrix

    """
    I = np.identity(3)
    ang = angle
    if degrees:
        ang = np.radians(ang)
    cos_id = np.cos(ang) * I
    print("cos_id",cos_id)
    sin_cross = np.sin(ang) * cross_product_matrix(axis_vector)
    print("sin_cross",sin_cross)
    cos_tens = (1-np.cos(ang)) * np.outer(axis_vector,axis_vector)
    print("cos_tens",cos_tens)
    rot_mat = cos_id + sin_cross + cos_tens

    print(rot_mat)

    return rot_mat
