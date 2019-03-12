Program overview and install
############################


Features
========

* Cross-program ONIOM calculation
* Interface with Gaussian, Turbomole, Molcas
* Reads output from CP2K, Quantum Espresso
* Mechanical, electrostatic, Ewald and Self-Consistent Ewald embedding
* Unique dimer detection
* Excitonic coupling via diabitazation
* Voronoi volume evaluation and visualisation

Requirements
============

* UNIX type system
* Python 2.7+ or 3.3+
* numpy (installed automatically)
* scipy (installed automatically)
* SWIG
* `Ewald <https://github.com/Crespo-Otero-group/Ewald>`_ (custom fork; only necessary for Ewald embedding)

Installation
============

1. Clone the repository to wherever you want to install it:

.. code-block:: bash

   cd /path/to/dir/
   git clone https://github.research.its.qmul.ac.uk/btx156/fromage.git
   cd fromage/

2. Install

.. code-block:: bash

  sudo pip install .

3. Set your environment variables

In your ``.bashrc``, add

.. code-block:: bash

  export FRO_GAUSS=g16
  export FRO_EWALD=Ewald

If you are using different binaries for Gaussian or Ewald, change accordingly.

Voilà!

