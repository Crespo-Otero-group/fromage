Miscellaneous scripts
#####################

Assign charges
==============

For the assignation of charges to moleculecular structures, cryspy uses their
connectivity. This is achieved through a set of tools all based on the
definition of a bond from just the information in an .xyz file. Two atoms are
considered bonded if their distance is below a certain threshold.

Armed with this approximated connectivity information, the charges assigned to
the atoms in one molecule can be redistributed to a whole cluster comprised of
the same molecule. One important approximation is that atoms with equivalent
connectivity (say the hydrogens in a methyl group) will have the same value: an
average of the three original ones.

This redistribution may be useful for other applications outside of cryspy, for
instance for the assignation of forcefield charges. To address this, the module
``assign_charges.py`` can be called as a script with as arguments the Gaussian
.log file containing the molecular charges, followed by a .xyz file containing
the cluster of molecules which need assignation.

In other words, the usage is:

.. code-block:: bash

  assign_charges.py mol.log clust.xyz

The output file will be called ``out_char`` by default. Several options are
available including specification of the maximum bond length or the type of
charge (ESP or Mulliken). More information can be found using:

.. code-block:: bash

  assign_charges.py --help

Population statistics
=====================

It is also often useful to quickly get some statistical information from a
population analysis done in Gaussian. For this, ``pop_stat.py`` can help.

Simply:

.. code-block:: bash

  pop_stat.py mol.log

This will output the maximum and minimum charge as well as the average absolute
chargem the standard deviation and the total calculated energy.

Pick a molecule
===============

This script selects a number of molecules from an .xyz file by indicating the
atom label of on of the constituent atoms. Something similar can often easily be
achieved by a visual program such as Chemcraft or VESTA, however this has the
advantage of working through a terminal if it becomes necessary:

.. code-block::

   pick_mol.py clust.xyz 7 13 42

Extract dimers
==============

TODO

Diabatization
=============

TODO
or VESTA
