Minimal Energy Conical Intersection optimisation
================================================

To optimise conical intersection geometries, the penalty function method of
Levine\ :cite:`Levine2008` is used, removing the need for nonadiabatic coupling
vectors. A function of the averaged S\ :sub:`1` and S\ :sub:`0` energies
(:math:`\bar{E}_{1-0}`) and the S\ :sub:`1`-S\ :sub:`0` energy gap
(:math:`\Delta E`) is minimised:

.. math::

  F = \bar{E}_{1-0} + \sigma \frac{\Delta E^2}{|\Delta E| + \alpha}

:math:`\sigma` is a Lagrangian multiplier and :math:`\alpha` is a parameter such
that :math:`\alpha \ll \Delta E`.


