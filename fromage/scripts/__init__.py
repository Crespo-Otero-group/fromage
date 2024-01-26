"""
Scripts

All of the executable files in fromage.

Modules
-------
    assign_charges
        Contains all of the routines for determining the connectivity of a group
        of atoms. Uses a gaussian.log file to 'populate' a .xyz file
    dimer_select
        Selects unique dimers from a .xyz file and writes them out
    pick_mol
        Extracts single moleculs from a .xyz file by one of their constituent
        atom's labels
    pop_stat
        General statistical analysis of a population analysis from Gaussian
    prepare_calculation
        Generates all of the necessary initial files for a fromage optimisation.
        This may include several calls to Ewald, which needs to be in the $PATH
    run_fromage
        Geometry optimisation via energy minimisation or penatly function of the
        gap

"""
#from fromage.scripts import *
