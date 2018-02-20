Program overview and install
############################


Features
========

* Cross-program ONIOM calculation [dapprich1999]_
* Interface with Gaussian, Trubomole, Molcas
* Reads output from CP2K, Quantum Espresso
* Mechanical, electrostatic, Ewald and Self-Consistent Ewald embedding
* Unique dimer detection
* Excitonic coupling via diabitazation
* Voronoi volume evaluation and visualisation

Requirements
============

* UNIX type system
* Python 2.7+ or 3.3+
* numpy
* scipy
* SWIG
* Ewald (custom fork)

Installation
============

1. Clone the repository to wherever you want to install it:

.. code-block:: bash

   cd /path/to/dir
   git clone https://github.research.its.qmul.ac.uk/btx156/cryspy.git

2. Compile

.. code-block:: bash

  cd cryspy/fdist
  swig -c++ -python fdist.i
  cd ../..
  python setup.py build_ext --inplace install

Voilà!

