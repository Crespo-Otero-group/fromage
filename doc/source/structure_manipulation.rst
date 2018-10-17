Structure manipulation in fromage
################################

As well as offering cluster model geometry optimisation, **fromage** can be used as a
library for the manipulation of molecular structures using Python. The main
vehicle for this is the Atom objects.

The Atom object
===============

The Atom object is a variable containing information about an atom in the
system, e.g. position, element or charge. In practice, one often deals with
molecules or unit cells as opposed to standalone atoms. Therefore many of the
functions in ``handle_atoms`` for instance, deal with lists of Atom objects in
their input.

IO
==

The functions for reading files can be found in ``read_file`` and the ones for
writing and editing in ``edit_file``. The most common file format for the
positions of atoms is .xyz which can be read by ``read_pos()`` for the final set
of coordinates or ``read_xyz()`` for all sets of coordinates in a list of lists.
``write_xyz()`` will write a list of atoms to an .xyz file.
