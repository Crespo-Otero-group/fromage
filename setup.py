#!/usr/bin/env python

from distutils.core import setup, Extension

fdist_module = Extension('_fdist', sources=['fdist_wrap.cxx', 'fdist.cpp'],)

setup(name='cryspy',
      version='0.9',
      author='Miguel Rivera, Michael Dommett',
      author_email='m.rivera@qmul.ac.uk',
      scripts=['assign_charges', 'cryspy', 'dimer_select',
               'pick_mol', 'pop_stat', 'prepare_calculation'],
      ext_modules=[fdist_module],
      py_modules=['assign_charges', 'atom', 'calc', 'cryspy', 'dimer_select', 'edit_file', 'fdist', 'handle_atoms',
                  'parse_config_file', 'periodic', 'pick_mol', 'pop_stat', 'prepare_calculation', 'read_file', 'volume'],
      )
