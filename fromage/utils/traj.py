""" The class for ordered successive steps of Mols

This can be used for MD trajectories, or optimisation steps.
"""
from fromage.utils.mol import Mol
import copy
import fromage.io.edit_file as ef
# the following looks like it is unused, but makes Traj list-like
import fromage.utils.listyness as listy # lgtm [py/import-and-import-from] lgtm [py/unused-import]
class Traj(object):
    """
    Class representing a list of Mols.

    This class uses composition to represent an ordered sequence of Mol
    objects, which themselves can be molecules, unit cells, etc.

    Attributes
    ----------
    frames : list of Mol objects
        The elements should be ordered
    """
    from fromage.utils.listyness import append, extend, insert, remove, index, pop, clear, count, __len__, __getitem__, __setitem__, __contains__

    def __init__(self, in_frames=[]):
        # if the user feeds one frame
        if isinstance(in_frames, Mol):
            in_frames = [in_frames]
        self.frames = in_frames
        self.chief_list = 'frames' # required to import listyness

    def __repr__(self):
        out_str = "Trajectory of " + str(len(self.frames)) + " frames."
        return out_str

    def __str__(self):
        return self.__repr__()

    def __add__(self, other_traj):
        new_traj = copy.deepcopy(self)
        new_traj.frames += other_traj.frames
        return new_traj

    def write_xyz(self, name):
        """Write an xyz file containing all geometries"""
        ef.write_traj_xyz(name, self)

