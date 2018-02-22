Input file description
######################

cryspy uses input files at two points of its execution. During the preparatory
calculation (using `prepare_calculation.py`), the `config` file is required.
During the actual geometry optimisation (`run_cryspy.py`), the file `cryspy.in`
is read if present.

config file
===========

All cryspy input files follow the same input structure. The order of the
keywords is irrelevant and blank lines are ignored. The keyword is stated and
then its value(s) after any number of whitespaces. Therefore:
.. code-block:: bash

  max_bl 1.9

Is the same as:

.. code-block:: bash

  max_bl      1.9

Below are listed the most important keywords available.

* `name`
  The name of your calculation

* `a_vec` , `b_vec` and `c_vec`
  The lattice vectors in Angstrom. For example:
.. code-block:: bash
  a_vec        8.9638004303         0.0000000000         0.0000000000
  b_vec        0.0000000000        10.5200004578         0.0000000000
  c_vec       -3.8748910079         0.0000000000        10.7924653741

* `cell_file`
  The file containing the atomic positions in the unit cell in .xyz format

* `high_pop_file`
  The file containing the population analysis used for the embedding of `mh`.
  Default: `gaussian_h.log`

* `high_pop_program`
  The program used to calculate the above file.Default: `gaussian`

* `high_pop_method`
  The method of population analysis. `Mulliken` or `ESP`. Default: `ESP`
