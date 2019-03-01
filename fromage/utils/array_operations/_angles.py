import numpy as np

def vec_angle(vector_1, vector_2, degrees = True):
    """
    Return the angle between two numpy vectors.

    The angle is brought into the range [-180,180] for degrees or [-1,1] for
    radians

    Parameters
    ----------
    vector_1, vector_2 : N x 1 numpy array
        The vectors whose angle needs to be calculated.
    degrees : bool (optional)
        Result in degrees or in radians. Default = True, so degrees
    Returns
    -------
    out_angle : float
        The angle between vectors
    """
    norm_1 = np.linalg.norm(vector_1)
    norm_2 = np.linalg.norm(vector_2)
    dot = np.dot(vector_1,vector_2)

    ang = np.arccos(dot/(norm_1*norm_2)) % (2 * np.pi)

    if degrees:
        ang = np.degrees(ang)
    return ang
