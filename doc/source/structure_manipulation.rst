Structure manipulation in fromage
#################################

As well as offering cluster model geometry optimisation, **fromage** can be used as a
library for the manipulation of molecular structures using Python.

To this end, importing fromage like e.g. ``import fromage as fro`` supplies the
user with a few useful functions such as ``fro.mol_from_file`` or
``fro.dimer_from_file`` which allow, in turn, grant the access to a host of
useful member functions. To understand how these functions operate, first it is
useful to get familiar with the data structure which is used to represent the
atoms of the system, the Atom object.

The Atom object
===============

The Atom object is a contains information about an atom in the system, e.g.
position, element or charge. In fact, the Atom object can also represent point
charges or any point in space, but as its name reflects, the bulk of its methods
are only useful when a chemical system is being represented. In practice, one
often deals with molecules or unit cells as opposed to standalone atoms. For
this, the Mol object should be employed.

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

The Mol object has an attribute ``mol.geom`` of the type ``GeomInfo``. This is
initialised as empty by default to save on computing effort, however useful
information can be requested such as the atomic positions in numpy array form
using ``mol.geom.coord_array()``. There are also functions to find the plane
which best contains the atomic coordinates with ``mol.geom.plane_coeffs()`` or three
orthonormal vectors giving a fingerprint to the molecular orientation:
``mol.geom.axes()``.

The Dimer object
================

Dimer objects represent pairs of molecules, each one being a Mol object stored
in the attributes ``dimer.mol_a`` and ``dimer.mol_b``. Since Mol objects can be
characterised by three axes, the angle between the axes of the two monomers can
supply a fingerprint for the dimer. They are supplied by the function
``dimer.angles()``.

IO
==

The functions for reading files can be found in ``read_file`` and the ones for
writing and editing in ``edit_file``. The most common file format for the
positions of atoms is .xyz which can be read by ``read_pos()`` for the final set
of coordinates or ``read_xyz()`` for all sets of coordinates in a list of lists.
``write_xyz()`` will write a list of atoms to an .xyz file.
