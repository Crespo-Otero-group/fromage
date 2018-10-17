"""
Fast distance

A distance calculator in C++. On certain machines, this can speed up the
calculation of a proximity volume by three. For optimal performance, use
fdist2 (the distance squared) rather than fdist if possible. This avoids the
evaluation of unnecesary square roots.

Module
------
    fdist
        Contains _fdist and _fdist2 for fast distance calculations
"""
#from cryspy.fdist import fdist
