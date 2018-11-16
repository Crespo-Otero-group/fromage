Ewald point charge fitting
==========================

Traditional electrostatic embedding for cluster models truncates the terms of
the Madelung summation up to the range of the cluster. This can sometimes be ill
advised since the Madelung summation is known to be slowly and conditionally
convergent.\ :cite:`Kittel1986`

A fast convergent reformulation of the Madelung summation is the Ewald
summation, expressed as a sum of both real and reciprocal space terms:

.. math::

  V^{Ewald} ( \mathbf{r} ) = \sum_{\mathbf{L}s} q_s \frac{\mathrm{erfc}{}(\gamma|\mathbf{r} - \mathbf{L}
  - \mathbf{R}_s|)} {|\mathbf{r} - \mathbf{L} - \mathbf{R}_s|}
  + \frac{4 \pi} {v_c}
  \sum_{\bf{G}\neq 0} \frac{1}{G^2} e^{-G^2/4\gamma^2}
  \Bigg{[}\sum_s q_s e^{i\mathbf{G}(\mathbf{r} - \mathbf{R}_s )}\Bigg{]}

Where :math:`\mathbf{L}` are the lattice translations, :math:`\mathrm{erfc}` is
the complementary error function, :math:`\gamma` is the arbitrary Ewald
constant, :math:`\mathbf{R}_s` the atomic sites in the unit cell, :math:`v_c`
the unit cell volume, :math:`\mathbf{G}` the reciprocal lattice translations and
:math:`q_s` the partial charge of the atoms in the unit cell.\
:cite:`Kantorovich2004`

To generate point charges fitted to an Ewald potential, **fromage** relies on
the program ``Ewald``, developed by Klintenberg, Derenzo and Weber.\
:cite:`Klintenberg2000,Derenzo2000`

In this program, an a supercell of point charges is generated. The Ewald
potential is computed in a central region and the points in the outer region are
fitted to reproduce said potential via direct summation. See citations above for
more details.

Ewald point charge embedding has successfully been used to describe excited
states in molecular crystals.\ :cite:`Dommett2017c,Wilbraham2016a,Presti2017`
