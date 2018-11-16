ONIOM cluster models
####################

Embedded Cluster
================

In order to investigate the local photochemical behaviour of molecular crystals,
it is practical to use a multiscale method which partitions the system in an
excited and a ground state region.

First, a molecule (or several) is selected from the unit cell to become the
:term:`model system<Model system>`. It is placed in the middle of a large
supercell from which a surrounding shell of molecules is selected, producing a
cluster (the :term:`real system<Real system>`). A molecule is selected if any of
its atoms fall within a chosen distance of the centroid of the model system.

With this partitioned cluster, we use the :term:`ONIOM` energy expression\
:cite:`Dapprich1999` to recover a multiscale energy:

.. math::

  E_{ONIOM} = E_{mh} + E_{rl} - E_{ml}

For straightforward :term:`mechanical embedding<Mechanical embedding>`, the
terms :math:`E_{mh}`, :math:`E_{rl}` and :math:`E_{ml}` are simply the energy of
the :term:`model system<Model system>` at a high level of theory, the
:term:`real system<Real system>` at a low level of theory, and the :term:`model
system<Model system>` at a low level of theory.

However this scheme only includes intersystem interactions at the low (ground
state) level of theory. To include intersystem Coulombic interactions in the
excited state, we use point charges. This is called :term:`electrostatic
embedding<Electrostatic embedding>` in the general :term:`ONIOM` literature and
the :term:`Embedded Cluster (EC) model<EC>` in this implementation.

The :term:`mh` calculation is embedded in point charges located at the atomic
sites of the surrounding cluster molecules. The value of these point charges is
not uniquely defined and a discussion is offered in a different
:ref:`section<pop>`.  In order to avoid the double counting of electrostatic
interactions, point charges must also be included in the :term:`ml` term.

Embedded cluster models have successfully been used by Ciofini's group in
characterising excited states in molecular crystals.\
:cite:`Presti2014,Presti2016,Presti2016a`

Ewald Embedded Cluster
======================

The :term:`EC` method described above represents short range interactions to a
reasonable extent. However the long range interactions, which are predominantly
electrostatic, are completely omitted. In fact they cannot be approximated by
increasing the size of the real system because the Madelung sum is conditionally
convergent\ :cite:`Kittel1986`. To remedy this, we use an Ewald embedding
scheme\ :cite:`Klintenberg2000,Derenzo2000` in :term:`mh` where a large array of
point charges at lattice positions is generated and then fitted to match the
Ewald potential. The combination of point charge Ewald embedding and :term:`EC`
is the :term:`Ewald Embedded Cluster (EEC) method<EEC>`.

This method requires some justification. First of all, the long range
electrostatic charges of the crystal are not cancelled in the :term:`ml` term.
If we wished, we could embed :term:`rl` and :term:`ml` in Ewald fitted point
charges. However when we perform geometry optimisation, the surrounding cluster
is fixed in place. Therefore the additional computation of the Ewald point
charges in the ground state Hamiltonians would only add a correcting constant
term.

Another first-glance objection is that the :term:`mh` charges from the
surrounding molecules which were included in the :term:`EC` model have
potentially been modified to match the Ewald potential, thus rendering the
cancellation of ground state interactions by the embedding of :term:`ml`
inexact. However by definition the Ewald potential contains the totality of the
Coulombic potential of the crystal, both short and long ranged. Furthermore a
spherical region of the Ewald point charge array of a chosen radius can be
chosen to remain of fixed charge, providing a 'buffer zone' from any highly
deviated charges which might break the point charge approximation.

Self-Consistent Ewald Embedded Cluster
======================================

A major omission from the :term:`EC` and :term:`EEC` models is the electrostatic
response to the excitation of the :term:`model system<Model system>` by the
surrounding cluster. To recover mutual polarisation effects, we employ an
extreme model where the entire crystal is excited at an electrostatic
equilibrium. The model system is embedded in Ewald point charges as in
:term:`EEC` at the optimised ground state position. A population analysis is
carried out on the model system whose charges are then redistributed in the
embedding supercell and gain fitted to the Ewald potential. This loop is
repeated until self-consistency.

The method is adapted from the work of Wilbraham *et al*.\
:cite:`Wilbraham2016a,Presti2017` Self-consistent Ewald embedding schemes were
previously used in the determination of NMR parameters\ :cite:`Weber2010`

This scheme is termed the Self-Consistent Ewald Embedded Cluster (\
:term:`SCEEC` S\ :sub:`1`). It accurately represents short range electrostatic
interactions from a mutually polarising delocalised excitation. Alternatively,
the self consistent loop can be performed in the ground state which would give
similar results to :term:`EEC` (:term:`SCEEC` S\ :sub:`0`).
