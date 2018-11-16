Interfaces
==========

For geometry optimisation, **fromage** communicates between different quantum
chemistry codes. They can be selected in the `fromage.in` by using the
`high_level` and `low_level` keywords. The current list of available programs
is:

 * Gaussian ``gaussian``\ :cite:`g09`
 * Molcas ``molcas``\ :cite:`Aquilante2016a`
 * Turbomole ``turbomole``\ :cite:`TURBOMOLE`
 * DFTB+ ``dftb``\ :cite:`Aradi2007` (under development)
