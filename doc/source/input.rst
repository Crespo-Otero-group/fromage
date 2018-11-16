Input file description
######################

**fromage** uses input files at two points of its execution. During the preparatory
calculation (using ``prepare_calculation.py``), the ``config`` file is required.
During the actual geometry optimisation (``run_fromage.py``), the file ``fromage.in``
is read if present.

config file
===========

All **fromage** input files follow the same input structure. The order of the
keywords is irrelevant and blank lines are ignored. The keyword is stated and
then its value(s) after any number of whitespaces. Therefore:

.. code-block:: python

  bond_thresh 1.9

Is the same as:

.. code-block:: python

  bond_thresh      1.9

Below are listed the most important keywords available.

name
  The name of your calculation. Default: ``fromage_calc``

a_vec, b_vec and c_vec
  The lattice vectors in Angstrom. This section has no default value. Example:

.. code-block:: python

  a_vec        8.9638004303         0.0000000000         0.0000000000
  b_vec        0.0000000000        10.5200004578         0.0000000000
  c_vec       -3.8748910079         0.0000000000        10.7924653741

vectors_file
  Alternatively the vectors can be stored in a file (called e.g. ``vectors``) of the form:

.. code-block:: python

  8.9638004303         0.0000000000         0.0000000000
  0.0000000000        10.5200004578         0.0000000000
  -3.8748910079         0.0000000000         10.7924653741

In which case the file name should be specified in the config file:

.. code-block:: python

  vectors_file vectors

cell_file
  The file containing the atomic positions in the unit cell in .xyz format.
  Default: ``cell.xyz``

high_pop_file
  The file containing the population analysis used for the embedding of ``mh``.
  Default: ``gaussian_h.log``

high_pop_program
  The program used to calculate the above file. Default: ``gaussian``

high_pop_method
  The method of population analysis. "Mulliken" or "ESP". Default: ``ESP``

low_pop_file
  The file containing the population analysis used for the embedding of ``ml``.
  Default: ``gaussian_l.log``

low_pop_program
  The program used to calculate the above file. Default: ``gaussian``

low_pop_method
  The method of population analysis. "Mulliken" or "ESP". Default: ``ESP``

bond_thresh
  The distance between two atoms in Angsrom below which **fromage** will consider the
  atoms to be bonded together. The definition of bond_thresh can be altered by
  using the keyword ``bonding``. Default: ``1.7``

bonding
  The method which determines whether two atoms are bonded. The options are
  ``dis``, ``cov`` and ``vdw``. ``dis`` measures the distance between two nuclei
  whereas ``cov`` measures the distance from the edge of the spheres of covalent
  radius and ``vdw`` from the edge of the sphere of van der Waals radius.
  Default: ``dis``

atom_label
  The number of the atom in the ``cell_file`` which belongs to the molecule which
  will become the :term:`model system<Model system>`. Several atoms can be
  specified and must be separated by whitespaces, however they must not belong
  to the same molecule. Default: ``1``

ewald
  Whether or not to use the Ewald embedding. To turn off, use "false", "no",
  "off", "zero", "none" or "nan" or any capitalisations. Default: ``off``

nchk
  The number of random points sampled around the model system by ``Ewald`` to
  check the accuracy of the fit. Default: ``1000``

nat
  The number of atoms included in the fixed charge region generated spherically
  by Ewald. Default: ``500``

an, bn and cn
  Multiplications of the unit cell along each cell direction to generate the
  Ewald supercell. The cell is multiplied 2N times per direction (N positive and
  N negative). Default: ``2``, ``2`` and ``2``

clust_rad
  The radius in Angstrom used to generate the cluster which will constitute the
  :term:`real system<Real system>`. The cluster includes all molecules which fit
  completely within the radius. Default: ``5``

self_consistent
  Whether or not to use the Self Consistent Ewald Embedding. Be sure to also
  turn on ``ewald``. Default: ``off``

sc_temp
  The template file for the self consistent population analyses. Default:
  ``sc_temp.template``

dev_tol
  The convergence threshold for the self consistent loop. This corresponds to
  the average deviation between two successive steps of the loop. Units in
  :math:`e^-` Default: ``0.001``

damping
  Damping factor for the self-consistent loop to solve certain convergence
  problems. Choose a value between 0 to 1 with 0 being no damping and 1 being
  complete damping (won't get you anywhere). Default: ``0``

print_tweak
  Whether or not to print the tweaked version of the cell with the selected
  molecule(s) completed and the whole cell centred around its centroid. This is
  useful for debugging and more involved analysis. Default: ``off``

fromage.in file
===============

The input structure is the same as for ``config``.

mol_file
  File name for the .xyz file containing the inital position of the :term:`model
  system<Model system>`. Default: ``mol.init.xyz``

shell_file
  File name for the .xyz file containing the molecules surrounding the
  :term:`model system<Model system>`. Default: ``shell.xyz``

out_file
  File name for the output file containing the geometry optimisation
  information. Default: ``fromage.out``

bool_ci
  Whether or not to optimise for :term:`MECI`. "1" for yes "0 for no. Default:
  ``0``

sigma
  The Lagrangian multiplier for the penalty function method for the location of
  :term:`MECI`. Only use if ``bool_ci`` is on. Defualt: ``3.5``

high_level
  The program used for the high level calculation. The options are ``gaussian``,
  ``dftb``, ``turbomole`` and ``molcas``. Default: ``gaussian``

low_level
  The program used for the low level calculation. The options are ``gaussian``,
  ``dftb``, ``turbomole`` and ``molcas``. Default: ``gaussian``

