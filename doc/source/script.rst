Miscellaneous scripts
#####################

Assign charges
==============

For the assignation of charges to molecular structures, fromage uses their
connectivity. This is achieved through a set of tools all based on the
definition of a bond from just the information in an .xyz file. Two atoms are
considered bonded if their distance is below a certain threshold.

Armed with this connectivity information, the charges assigned to the atoms in
one molecule can be redistributed to a whole cluster comprised of the same
molecule. One important approximation is that atoms with equivalent connectivity
(say the hydrogens in a methyl group) will have the same value: an average of
the three original ones.

This redistribution may be useful for other applications outside of fromage, for
instance for the assignation of forcefield charges. To address this, the module
``fro_assign_charges.py`` can be called as a script with as arguments the
Gaussian .log file containing the molecular charges, followed by a .xyz file
containing the cluster of molecules which need assignation.

In other words, the usage is:

.. code-block:: bash

  fro_assign_charges.py mol.log clust.xyz

The output file will be called ``out_char`` by default. Several options are
available including specification of the maximum bond length or the type of
charge (ESP or Mulliken). More information can be found using:

.. code-block:: bash

  fro_assign_charges.py --help

Population statistics
=====================

It is also often useful to quickly get some statistical information from a
population analysis done in Gaussian. For this, ``fro_pop_stat.py`` can help.

Simply:

.. code-block:: bash

  fro_pop_stat.py mol.log

This will output the maximum and minimum charge as well as the average absolute
charge, the standard deviation and the total calculated energy.

Pick a molecule
===============

This script selects a number of molecules from an .xyz file by indicating the
atom label of on of the constituent atoms. Something similar can often easily be
achieved by a visual program such as Chemcraft or VESTA, however this has the
advantage of working through a terminal if it becomes necessary:

.. code-block:: bash

   fro_pick_mol.py clust.xyz 7 13 42

This will pick out the molecules containing atoms 7, 13 and 42 from the cluster
of molecules ``clust.xyz``.

Manipulate unit cell .xyz files
===============================

``fro_uc_tools.py`` can be used to produce supercells, molecular clusters and other
ouputs based on the unit cell. It requires two inputs: a unit cell file in .xyz
format and a vectors file like this:

.. code-block:: python

  8.9638004303         0.0000000000         0.0000000000
  0.0000000000        10.5200004578         0.0000000000
  -3.8748910079         0.0000000000         10.7924653741

There are many options to choose from so it is suggested to use ``fro_uc_tools -h``
for instructions.

For example:

.. code-block:: bash

    fro_uc_tools.py cell.xyz vectors -r 12

Will use a unit cell geometry file (``cell.xyz``) and its associated lattice
vectors (``vectors``) to produce a molecular cluster of all complete molecules
with atoms falling within a radius of 12 Angstroms from the origin (``-r 12``).
The output is called ``cluster_out.xyz`` by default.

Analyse dimers in aggregate geometries
======================================

Whether one is presented with a cell stemming from a periodic boundary condition
calculation, or an oligomer in a localised cluster, the geometric properties of
the dimers present within can help elucidate some of the intermolecular features
of the system. The script ``fro_dimer_tools.py`` can identify the unique dimers
in the supplied geometry, taking into account periodicity if relevant. The
dimers can further be characterised by the angles between their principal,
secondary and normal axes, as well as their centroid-to-centroid distance. As
before, many parameters can be altered so using ``fro_dimer_tools.py -h`` is
encouraged.

A suggested use is the following:

.. code-block:: bash

    fro_dimer_tools.py clust.xyz -v -p

Will analyse an ``.xyz`` geometry file of a cluster of molecules
(``clust.xyz``), with a verbose output (``-v``), and print (``-p``) all of the
unique dimers it finds within the cluster. A file called ``dimers.dat`` will
also be printed with some geometric information related to each dimer and a
suggested classification, being ``S-S``, ``E-F`` or ``F-F`` (side-by-side,
edge-to-face or face-to-face).

Voronoi volume evaluation
=========================

It can be useful to determine the available volume of a molecule in an aggregate
environment. To do this, one could use the union of the van der Waals volumes of
each atom, or the Voronoi volume of the molecule, scaled by van der Waals radii.

.. code-block:: bash

    fro_volumetrics.py clust.xyz -l 13

This will produce cube files of the available volume of the molecule containing
atom 13 (``-l 13``) within the cluster of molecules (``clust.xyz``). The ouputs
are the Voronoi volume (``voro.cube``), the van der Waals volume (``vdw.cube``)
and the union of the two (``add.cube``). A file called ``volumes`` prints the
integrated volume of each of the three.

Exciton coupling evaluation
===========================

Exciton coupling evaluation from Gaussian output files can also be carried out,
using ``fro_coupling.py``. A diabatisation of the Hamiltonian is employed which
relies on the calculation of excited state properties such as population
analysis or transition dipole moments.\ :cite:`Arag2015` More options are also
available.

As an example of use, the line:

.. code-block:: bash

    fro_coupling.py -m DIA -p TDM -mf a.log b.log -of dim_ab.log -os 2

Will use the diabatisation method (``-m DIA``) and use the transition dipole
moment property (``-p TDM``) to read the Gaussian log files of monomer S\
:sub:`1` calculations (``-mf a.log b.log``) and the dimer S\ :sub:`2`
calculation (``-of dim_ab.log``) with state of interest S\ :sub:`2` (``-os 2``).
The output will show the diabatic Hamiltonian, whose off-diagonal elements are
the exciton coupling values.

Exciton classification
======================

It is sometimes useful to classify excitons as localised, delocalised or
charge-transfer. To this end, the script ``fro_exciton_classification.py`` uses
a Mulliken partition scheme to analyse the migration of charge density within a
single excitation from a TDDFT or CIS calculation.\ :cite:`Crespo-Otero2012` :cite:`Sen2013`

To use this, prepare an Gaussian calculation of a dimer in a given excited state
using either TDDFT or CIS. Make sure that all of the atoms of one molecule (A)
appear before all of the atoms of the second molecule (B). Also make sure to
print the ``.rwf`` file by using the input tag ``%rwf=filename.rwf``.

Here is an example Gaussian input for a range-separated hybrid TDDFT calculation
up to the fourth excited state:

.. code-block:: html

    %chk=title.chk
    %mem=[X]GB
    %nproc=[X]
    #p wb97xd 6-31g* gfprint pop=full

    title

    0 1
         H   0.00 0.00 0.00
         C   0.00 0.00 0.00
         .
         .
         .

    --link1--
    %chk=title.chk
    %rwf=title
    #p wb97xd 6-31g* td(nstates=4,root=4) gfprint pop=full IOP(3/33=3) geom=allcheck guess=read density=current

Now run:

.. code-block:: bash

    fro_exciton_classification.py tddft.log tddft.rwf 3

This reads first the Gaussian log file (``.log``), then the read-write file
(``rwf``) and analysises the third excited state (``3``). Two indices relating
to the electron density migration will be printed, and a classification of the
exciton will be suggested as ``LOC(A)``, ``LOC(B)``, ``CT A->B``, ``CT B->A`` or
``Delocalised`` (localised on A, or on B, charge transfer from A to B, vice
versa, or delocalised). More details can be found in the references above.

