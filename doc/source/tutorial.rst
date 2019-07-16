Tutorial for geometry optimisation
##################################

This tutorial will guide you step by step through a typical ONIOM calculation
done with **fromage**. This page is meant to be self contained but for a
refresher on the abbreviations used, visit the :ref:`glossary<gloss>`.
Everything discussed herein is elaborated upon in the rest of the documentation.
The example calculation is rather computationally costly and time intensive.
This is because we choose an model system with a distinctive excited state minimum
and conical intersection geometry.

The example for this tutorial is the investigation of the emission properties of
the (2-hydroxyphenyl)propenone (HP) crystal using the Ewald Embedded Cluster
model (EEC) and Gaussian. The extension of this protocol to other programs is
straightforward. The following steps will involve calculating the geometries and
associated energies of the ground state, first excited state and S\ :sub:`1` -
S\ :sub:`0` Minimal Energy Conical Intersection (MECI).

The necessary input files are in the ``tutorial`` directory. We follow the ONIOM
nomenclature for the different regions of the embedded cluster.  In other words
the whole cluster is the :term:`"real system"<Real system>` and the central
molecule is called the \ :term:`"model system"<Model system>`\ . We abbreviate
the principal calculations involved in the ONIOM method as ``mh`` for the
:term:`model system<Model system>` at a high level of theory, ``ml`` for the
:term:`model system<Model system>` at a low level of theory and ``rl`` for the
:term:`real system<Real system>` at a low level of theory. For the location of
conical intersections, a fourth calculation is typically involved, ``mg``, which
indicates the :term:`model system<Model system>` at the high level of theory but
in the ground state.

Calculation set up
==================

The initial step is to generate the Ewald charge background as well as all
of the necessary files needed for geometry optimisation.

Input
-----

Make a directory where the preparatory calculation will take place:

.. code-block:: bash

  mkdir opt_1/
  cd opt_1/
  cp path/to/fromage/tutorial/* .

The input files that you have just copied are:

* ``cell.xyz``
    A file containing the optimised positions of the atoms in
    the unit cell

* ``config``
    The **fromage** configuration file, described below

* ``high_pop.log``
    A Gaussian output file of a population analysis for one
    molecule using the high level of theory (in this case PBE 6-31G*)

* ``low_pop.log``
    The same but for the low level of theory (here, HF
    STO-3G)

* ``*.template``
    The template files for your upcoming ONIOM calculation.
    This is includes ``mh.template``, ``ml.template``, ``rl.template``
    and ``mg.template``. Note that the checkpoint file name is ``gck.chk``
    in all cases and that the levels of theory match those of the population log
    files

Execution
---------

If your installation was successful, all of the **fromage** scripts should be in your
system path already. In that case, running the program simply involves typing:

.. code-block:: bash

  fro_prep_calc.py

Output
------

After a few minutes, you will be greeted with a series of outputs:

* ``prep.out``
    Output file with some information about the setup
    calculation

* ``mol.init.xyz``
    The initial position of the :term:`model system<Model system>`

* ``shell.xyz``
    The molecules surrounding the :term:`model system<Model system>`

* ``mh/ ml/ rl/ mg/``
    Directories containing a ``.temp`` file each. For
    example ``mh/`` contains ``mh.temp``

* ``ewald/``
    The directory where the ewald calculation is run. The outputs in here are
    not important for this tutorial

Geometry optimisation
=====================

We will calculate the geometries and associated energies of the ground state
minimum, first excited state minimum and S\ :sub:`1` - S\ :sub:`0` Minimal
Energy Conical Intersection (MECI).

Ground state
------------

Input
^^^^^

These are all the files needed for the geometry optimisation. Most of them
were already generated from the previous step.

* ``fromage.in``
    The input file which contains the specifications for the geometry
    optimisation

* ``mol.init.xyz``
    See above

* ``shell.xyz``
    See above

* ``mh/ ml/ rl/``
    Directories containing their respective ``.temp`` files


Execution
^^^^^^^^^

An important part of calculations in **fromage** is the assignment of memory to each
component calculation. Some times, depending on the system size and the
combination of methods used, ``rl`` will need more memory than ``mh``. Make sure
to adapt the memory requested in all three ``.temp`` files to match the capacity
of your system.

When this is ready, submit your job with the command:

.. code-block:: bash

  fro_run.py

On the command line or in your job queue.

Output
^^^^^^

You can expect this calculation to take a few hours depending on your
computational resources. The convergence criterion of the optimisation is very
strict by default so it is up to the user's judgement whether they wish to abort
the calculation once they have achieved a satisfactory precision.

* ``fromage.out``
    The main output file. This contains information about the energies and
    gradients at each step of the optimisation

* ``geom_mol.xyz``
    The positions of the :term:`model system<Model system>` throughout the optimisation

* ``geom_clust.xyz``
    The position of the real system throughout the optimisation. Only the
    :term:`model system<Model system>` will change

``geom_mol.xyz`` should show very slight rearrangement of the molecule since its
Gaussian-optimised ground state geometry is close its crystal.

Vertical excitation
^^^^^^^^^^^^^^^^^^^

To calculate the vertical excitation, make a new directory called ``exci/`` and
copy the ``mh.com`` file to it (it should now contain the optimised geometry):

.. code-block:: bash

  mkdir exci/
  cp mh/mh.com exci/
  cd exci/

Then edit the ``mh.com`` file to remove the ``force`` keyword and add
``td(nstates=5,root=1)``.
Now run Gaussian (or use a submission script):

.. code-block:: bash

  g16 mh.com

This will give the excitation energies of the first
five excited states, easily accessible with a judicious ``grep``:

.. code-block:: bash

  grep 'Excited State' mh.log

The output will look like this:

.. code-block:: bash

  Excited State   1:      Singlet-?Sym    4.1338 eV  299.93 nm  f=0.6373  <S**2>=0.000
  Excited State   2:      Singlet-?Sym    4.3687 eV  283.80 nm  f=0.0068  <S**2>=0.000
  Excited State   3:      Singlet-?Sym    4.5466 eV  272.70 nm  f=0.0760  <S**2>=0.000
  Excited State   4:      Singlet-?Sym    5.4041 eV  229.43 nm  f=0.0318  <S**2>=0.000
  Excited State   5:      Singlet-?Sym    5.8322 eV  212.59 nm  f=0.0431  <S**2>=0.000

First excited state
-------------------

We now wish to optimise the geometry of the molecule in the first excited state.
The procedure will be almost identical to the one for the ground state
optimisation but with an added keyword to ``mh.temp``.

Input
^^^^^

First, copy the whole ``opt_1`` directory to conserve the ground state data.
Presuming you are still in ``opt_1/exci/``, just type:

.. code-block:: bash

  cd ../../
  cp -r opt_1/ opt_2/
  cd opt_2/

Now edit the file ``mh/mh.temp`` to add the keyword ``td(nstates=1,root=1)``.

And edit your ``mol.init.xyz`` to match the last geometry in ``geom_mol.xyz``
from your ``opt_1/`` directory.

Execution
^^^^^^^^^

As usual, type:

.. code-block:: bash

  fro_run.py

Or submit it to your job scheduler.

Output
^^^^^^

This should typically take longer than your ground state calculation if you have
succeeded in setting it up in a way that the limiting calculation is ``mh``.

As described above, you will receive ``fromage.out``, ``geom_mol.xyz`` and
``geom_clust.xyz``.

This time, you should be able to see the excited state proton transfer in
``geom_mol.xyz`` as the optimised structure is in keto form.

MECI
----

Input
^^^^^

One final time, copy the whole directory:

.. code-block:: bash

  cd ..
  cp -r opt_2/ opt_3/
  cd opt_3/

Edit the ``mol.init.xyz`` file with the final geometry of
``opt_2/geom_mol.xyz``.

And in ``fromage.in``, add a line at the bottom ``bool_ci 1``. This turns on MECI
search. Keep in mind that this calculation will use ``mg`` so change the memory
requested in all of your ``.temp`` files accordingly.

Execution
^^^^^^^^^

Again:

.. code-block:: bash

  fro_run.py

And wait however long it takes.

Output
^^^^^^

The usual ``fromage.out``, ``geom_mol.xyz`` and ``geom_clust.xyz`` will be
generated.

``fromage.out`` will contain different information, pertaining to the value and
gradients of the penalty function which is being minimised instead of the
energy.




