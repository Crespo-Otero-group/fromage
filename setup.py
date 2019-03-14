#!/usr/bin/env python

from setuptools import setup, Extension

fdist_module = Extension('fromage.fdist._fdist', sources=['fromage/fdist/fdist.i'])

setup(name='fromage',
      version='1.0',
      description='FRamewOrk for Molecular Aggregate Excitations',
      author='Miguel Rivera, Michael Dommett, Rachel Crespo-Otero',
      author_email='r.crespo-otero@qmul.ac.uk',
      license='MIT',
      ext_modules=[fdist_module],
      packages=['fromage',
                'fromage.fdist',
                'fromage.io',
                'fromage.scripts',
                'fromage.utils',
                'fromage.utils.exci_coupling',
                'fromage.utils.array_operations',
                'fromage.utils.mol'],
      scripts=['fromage/scripts/fro_assign_charges.py',
               'fromage/scripts/fro_run.py',
               'fromage/scripts/fro_dimer_tools.py',
               'fromage/scripts/fro_pick_mol.py',
               'fromage/scripts/fro_pop_stat.py',
               'fromage/scripts/fro_prep_run.py',
               'fromage/scripts/fro_coupling.py',
               'fromage/scripts/fro_volumetrics.py',
               'fromage/scripts/fro_uc_tools.py'],
      install_requires=[
          'numpy',
          'scipy',],
      )
