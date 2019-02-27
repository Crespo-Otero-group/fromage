import numpy as np

def coord_array(self):
    """
    Return a numpy array of the coordinates

    Returns
    -------
    arr : Nat x 3 numpy array
        Array of the form [[x1,y1,z1],[x2,y2,z2],...]

    """
    nat = len(self)
    arr = np.zeros((nat,3))
    for i,atom in enumerate(self):
        arr[i][0] = atom.x
        arr[i][1] = atom.y
        arr[i][2] = atom.z
    return arr

def pairwise_distances2(self):
    """
    Return a numpy array of the pairwise squared distances

    Returns
    -------
    arr : Nat x Nat numpy array
        Squared distance matrix (lower triangular and 0 diagonal)

    """

    nat = len(self)
    dis2 = np.zeros((nat,nat))
    arr = self.coord_array()
    for i, row_i in enumerate(arr):
        for j, row_j in enumerate(arr[i+1:]):
            squared_dist = np.sum(np.square(row_i-row_j))
            dis2[j+i+1][i] = np.sqrt(squared_dist)
            #print(np.sqrt(squared_dist))
    print(dis2)
    return
