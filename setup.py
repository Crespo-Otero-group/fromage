#!/usr/bin/env python

from distutils.core import setup, Extension

fdist_module = Extension('cryspy.fdist._fdist', sources=['cryspy/fdist/fdist_wrap.cxx', 'cryspy/fdist/fdist.cpp'],)

setup(name='cryspy',
      version='1.0',
      author='Miguel Rivera, Michael Dommett, Rachel Crespo-Otero',
      author_email='r.crespo-otero@qmul.ac.uk',
      ext_modules=[fdist_module],
      packages=['cryspy',
                'cryspy.fdist',
                'cryspy.io',
                'cryspy.scripts',
                'cryspy.utils'],
      scripts=['cryspy.scripts.assign_charges',
               'cryspy.scripts.cryspy',
               'cryspy.scripts.dimer_select',
               'cryspy.scripts.pick_mol',
               'cryspy.scripts.pop_stat',
               'cryspy.scripts.prepare_calculation'],
      )
