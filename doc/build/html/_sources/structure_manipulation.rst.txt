Structure manipulation in fromage
#################################

As well as offering cluster model geometry optimisation, **fromage** can be used as a
library for the manipulation of molecular structures using Python. The main
vehicle for this is the Atom objects.

The Atom object
===============

The Atom object is a contains information about an atom in the system, e.g.
position, element or charge. In practice, one often deals with molecules or unit
cells as opposed to standalone atoms. For this, the Mol object should be
employed.

The Mol object
==============

The Mol object is related to usual Python lists by composition. It is meant to
contain lists of Atom objects which represents points in space with associated
properties. The Mol object has a multitude of methods which are ready to use for
quick geometry manipulation.

Many of these methods take advantage of the modularity of molecular crystals, as
such the intermolecular distance which constitutes a bond is a crucial
parameter. The default parameters for this should usually be sufficient but in
case the geometry being investigated is particularly distorted, the bonding can
be defined in different ways: the distance between nuclei, between vdW radii and
covalent radii. ``set_bonding()`` is used to tweak this definition.

IO
==

The functions for reading files can be found in ``read_file`` and the ones for
writing and editing in ``edit_file``. The most common file format for the
positions of atoms is .xyz which can be read by ``read_pos()`` for the final set
of coordinates or ``read_xyz()`` for all sets of coordinates in a list of lists.
``write_xyz()`` will write a list of atoms to an .xyz file.
