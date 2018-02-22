.. _gloss:

Essential glossary
#########################

Throughout this documentation, we employ terms and abbreviations which may not
be immediately clear to the new user. Herein is a compiled list of important
terms accompanied by brief descriptions.  Furthur discussion for many of these
terms is included in the user documentation.


Cluster models
==============

In order to take into account the environmental energetical contributions to an
excited molecule inside a crystal, we offer several kinds of non-intrusive
embedding methods:

.. glossary::

  ONIOM
    Our own N-layered Integrated mlecular Orbital and Molecular
    mechanics. An extrapolative (subtractive) embedding paradigm for multilevel
    calculations

  EC
    Embedded Cluster. The application of the ONIOM method with
    electrostatic embedding to a cluster of molecules taken from their crystalline
    positions

  EEC
    Ewald Embedded Cluster. Similar to EC but using an electrostatic
    embedding scheme aimed at reproducing the Ewald potential of the crystal

  SCEEC
    Self-Consistent Ewald Embedded Scheme. The EEC model where the
    embedding charges are computed self consistently in the excited state

  Mechanical embedding
    ONIOM embedding where no point charges are included and the intersystem
    electrostatic interaction is purely at low level

  Electrostatic embedding
    ONIOM embedding where the point charges from the surrounding cluster are
    included in the Hamiltonian of the high level of theory to be treat
    Coulombic terms in the excited state

Calculations
============

In the expression of the 2-level ONIOM energy, the system under scrutiny is
partitioned into two regions:

.. glossary::

  Model system
    The central molecule(s) whose properties are under scrutiny

  Real system
    The complete cluster of molecules taken from their
    crystalline positions, including and surrounding the model system

By nature, 2-level ONIOM combines two levels of theory:

.. glossary::

  High level
    A level of theory which describes the required properties of
    the model system. In photochemistry this will often be an excited state method

  Low level
    A less computationally expensive level of theory than the high
    level. This will typically be a ground state method

When optimising the geometry of a molecule using one of the cluster models, a
few concurrent calculations must be carried out.

.. glossary::

  mh
    Model system high level calculation

  ml
    Model system low level calculation

  rl
    Real system low level calculation

  mg
    Model system calculation in the ground state of the high level

General photochemistry
======================

.. glossary::

  MECI
  Minimum Energy Conical intersection
