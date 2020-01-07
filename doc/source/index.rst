.. fromage documentation master file, created by
   sphinx-quickstart on Tue Feb 13 10:35:39 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: logo.png
   :width: 400 px

|

**fromage** (FRamewOrk for Molecular AGgregate Excitations) is a library designed
to support investigation of photophenomena in molecular crystals. Among the
features are:

 * Cross-program ONIOM-style calculations with different electrostatic embedding
   methods
 * Location of energy minima and minimal energy conical intersections
 * Command line geometry manipulation tools
 * Evaluation of exciton coupling values using multiple schemes (under
   development)

The current version is 1.0

To cite the use of the program, please use:

Rivera, M., Dommett, M., Sidat, A., Rahim, W., Crespo‐Otero, R. fromage: A library for the study of molecular crystal excited states at the aggregate scale. *J Comput Chem* 2020; 1– 14. https://doi.org/10.1002/jcc.26144

And if you are using one of the ONIOM implementations:

Rivera, M., Dommett, M., Crespo-Otero, R. ONIOM(QM:QM′) Electrostatic Embedding Schemes for Photochemistry in Molecular Crystals. *J. Chem. Theory Comput.* 2019; 15, 4, 2504-2516 https://doi.org/10.1021/acs.jctc.8b01180

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: General

   overview
   gloss
   tutorial
   license
   contact
   zrefs

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Theoretical background

   cluster_models
   ewald
   penalty_func
   pop

.. toctree::
  :maxdepth: 2
  :numbered:
  :caption: Program documentation

  interfaces
  input
  structure_manipulation
  script
  modules


References and indices
======================

* :ref:`refs`
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
