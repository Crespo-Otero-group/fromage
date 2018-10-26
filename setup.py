#!/usr/bin/env python

from distutils.core import setup, Extension

fdist_module = Extension('fromage.fdist._fdist', sources=['fromage/fdist/fdist_wrap.cxx', 'fromage/fdist/fdist.cpp'],)

setup(name='fromage',
      version='1.0',
      author='Miguel Rivera, Michael Dommett, Rachel Crespo-Otero',
      author_email='r.crespo-otero@qmul.ac.uk',
      ext_modules=[fdist_module],
      packages=['fromage',
                'fromage.fdist',
                'fromage.io',
                'fromage.scripts',
                'fromage.utils'],
      scripts=['fromage/scripts/assign_charges.py',
               'fromage/scripts/run_fromage.py',
               'fromage/scripts/dimer_select.py',
               'fromage/scripts/pick_mol.py',
               'fromage/scripts/pop_stat.py',
               'fromage/scripts/prepare_calculation.py',
               'fromage/scripts/uc_tools.py'],
      install_requires=[
          'numpy',
          'scipy',],
      )
