Tutorial
########

This tutorial will guide you step by step through a typical calculation done
with cryspy. This page is meant to function as a standalone but for a refresher
on the abbreviations used, visit the :ref:`glossary<gloss>`. Everything discussed
herein is elaborated upon in the rest of the documentation.

The model calculation used here is the investigation of the emission properties
of the (2-hydroxyphenyl)propenone (HP) crystal using the Ewald Embedded Cluster
model (EEC) and Gaussian.  This will involve calculating the geometries and
associated energies of the ground state, first excited state and S:sub:`1` -
S:sub:`0` Minimal Energy Conical Intersection (MECI).


The necessary input files are in the ``tutorial`` directory. We
follow the ONIOM nomenclature for the different regions of the embedded cluster.
In other words the whole cluster is the "real" system while the central
molecule is called the "model" system. We abbreviate the principal
calculations involved in the ONIOM method as ``mh`` for the model system at a
high level of theory, ``ml`` for the model system at a low level of theory
and ``rl`` for the real system at a low level of theory. For the location of
conical intersections, a fourth calculation is typically involved, ``mg``,
which indicates the real system at the high level of theory but in the ground
state.

To follow along, start from the ``tutorial`` directory.

Calculation set up
==================

The initial step step involves generating the Ewald charge background as well as
all of the necessary files needed for geometry optimisation.

Input
-----

Make a directory where the preparatory calculation will take place:

.. code-block:: bash

  mkdir prep
  cp * prep
  cd prep

The input files that you have just copied are:

* ``cell.xyz``
    A file containting the optimised positions of the atoms in
    the unit cell

* ``config``
    The cryspy configuration file, described below

* ``high_pop.log``
    A Gaussian output file of a population analysis for one
    molecule using the high level of theory (in this case PBE 6-31G*)

* ``low_pop.log``
    The same but for the low level of theory (here, HF
    STO-3G)

* ``*.template``
    The template files for your upcoming ONIOM calculation.
    this is includes ``mh.template``, ``ml.template``, ``rl.template``
    and ``mg.template``. Note that the checkpoint file name is ``gck.chk``
    in all cases and that the levels of theory match those of the population log
    files

Execution
---------

If your installation was successful, all of the cryspy scripts should be in your
system path already. In that case, running the program simply involves typing:

.. code-block:: bash

  prepare_calculation.py

Output
------

After a few minutes, you will be greeted with a series of outputs:

* ``prep.out``
    Output file with some information about the setup
    calculation

* ``mol.init.xyz``
    The initial position of the model system

* ``shell.xyz``
    The molecules surrounding the model system

* ``mh ml rl mg``
    Directories containing a ``.temp`` file each. For
    example ``mh`` contains ``mh.temp``

* ``ewald``
    The directory where the ewald calculation is run. The outputs in here are
    not important for this tutorial

Geometry optimisation
=====================

We will calculate the geometries and associated energies of the ground state,
first excited state and S:sub:`1` - S:sub:`0` Minimal Energy Conical
Intersection (MECI).

Ground state
------------

Input
^^^^^

These are all the files needed for the geometry optimisation. Most of them
were already generated from the previous step.

* ``cryspy.in``
    The input file which contains the specifications for the geometry
    optimisation

* ``mol.init.xyz``

* ``shell.xyz``

* ``mh ml rl`` directories containting their respective ``.temp`` files

Execution
^^^^^^^^^

An important part of calculations in cryspy is the assignement of memory to each
component calculation. Some times, depending on the system size and the
combination of methods used, ``rl`` will need more memory than ``mh``. Make sure
to adapt the memory requested in all three ``.temp`` files to match the capacity
of your system.

When this is ready, submit your job with the command:

.. code-block:: bash

  run_cryspypy

Output
^^^^^^

