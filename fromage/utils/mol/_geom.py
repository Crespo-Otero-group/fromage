import numpy as np
import fromage.utils.array_operations as ao

def coord_array(self):
    """
    Return a numpy array of the coordinates

    Returns
    -------
    coord_arr : Nat x 3 numpy array
        Array of the form [[x1,y1,z1],[x2,y2,z2],...]

    """
    nat = len(self)
    coord_arr = np.zeros((nat,3))
    for i,atom in enumerate(self):
        coord_arr[i][0] = atom.x
        coord_arr[i][1] = atom.y
        coord_arr[i][2] = atom.z
    return coord_arr
