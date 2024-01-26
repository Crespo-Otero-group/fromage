""" Opens xyz files and returns the atoms and cooridnates in various forms
"""
import numpy as np
def open_xyz(xyz_file):
    """
    Opens an xyz file and returns a list of lists

    Parameters
    ----------
    xyz_file: String
        File path

    Returns
    ----------
    List of lists: List of lists
        Outer list is whole xyz file, each inner list is a line of the file containing
        the symbol and x,y,z coordinates

    """
    file=open(xyz_file,'r').read().splitlines()
    # Check if first line contains number of atoms
    if len(file[0].split())!= 1:
        exit("Not a proper xyz file!")
    natoms = int(file[0])
    if natoms != len(file[2:]):
        exit("Number of atoms in header does not match number of atoms in file!")
    else:
        return file[2:]

def xyz_to_matrix(xyz_list):
    """
    Takes a list of lists containing xyz file lines and returns a coordinate matrix

    Parameters
    ----------
    xyz_list: List of lists
        Outer list is whole xyz file, each inner list is a line of the file containing
        the symbol and x,y,z coordinates

    Returns
    ----------
    coordinate_matrix: np.ndarray
        Nx3 matrix where N=number of atoms
    """
    coordinate_matrix=np.zeros((len(xyz_list),3))
    for i in range(len(xyz_list)):

        coordinate_matrix[i,0]=xyz_list[i][1]
        coordinate_matrix[i,1]=xyz_list[i][2]
        coordinate_matrix[i,2]=xyz_list[i][3]
    return coordinate_matrix

def symbols_from_xyz(xyz_list):
        """
        Takes a list of lists containing xyz file lines and the elemental symbols

        Parameters
        ----------
        xyz_list: List of lists
            Outer list is whole xyz file, each inner list is a line of the file containing
            the symbol and x,y,z coordinates

        Returns
        ----------
        symbols: List of strings
            List of atomic symbols
        """
        symbols=[]
        for i in range(len(xyz_list)):
            symbols.append(str(xyz_list[i][0]))
        return symbols
