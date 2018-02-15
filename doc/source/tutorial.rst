Tutorial
########

This tutorial will guide you step by step through a typical calculation done
with cryspy. This page is meant to function as a standalone but for a refresher
on the abbreviations used, visit the :ref:`glossary`. Everything discussed
herein is elaborated upon in the rest of the documentation.

The model calculation used here is the investigation of the emission properties of
the (2-hydroxyphenyl)propenone (HP) crystal using the Ewald Embedded Cluster
model (EEC) and Gaussian.

The necessary input files are in the :code:`tutorial` directory. We
follow the ONIOM nomenclature for the different regions of the embedded cluster.
In other words the whole cluster is the "real" system while the central
molecule is called the "model" system. We abbreviate the principal
calculations involved in the ONIOM method as :code:`mh` for the model system at a
high level of theory, :code:`ml` for the model system at a low level of theory
and :code:`rl` for the real system at a low level of theory. For the location of
conical intersections, a fourth calculation is typically involved, :code:`mg`,
which indicates the real system at the high level of theory but in the ground
state.

Calculation set up
==================

Necessary input files
---------------------

* :code:`cell.xyz`: A file containting the optimised positions of the atoms in
  the unit cell

* :code:`config`: The cryspy configuration file, described below

* :code:`high_pop.log`: A Gaussian output file of a population analysis for one
  molecule using the high level of theory (in this case PBE 6-31G*)

* :code:`low_pop.log`: The same but for the low level of theory (here, HF
  STO-3G)

* :code:`*.template`: The template files for your upcoming ONIOM calculation






