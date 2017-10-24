#!/usr/bin/env python
"""Reads the user's config inputs

This contains all of the default settings and may raise warnings and exceptions
if contradictory features are included.
"""
import read_file as rf

def parse_config(name = "config"):
    # Default parameters
    inputs = {
    "cell_file":"cell.xyz",
    "high_pop_program":"gaussian",
    "high_gauss_file":"gaussian.log",
    "high_gauss_method":"ESP",
    "high_cp2k_file":"cp2k.out",
    "high_cp2k_method":"ESP",
    "low_pop_program":"gaussian",
    "low_gauss_file":"gaussian.log",
    "low_gauss_method":"ESP",
    "low_cp2k_file":"cp2k.out",
    "low_cp2k_method":"ESP",
    "max_bl":"1.7",
    "atom_label":"1",
    "ewald":"1", # gets cast to int then bool
    "nChk":"1000",
    "nAt":"500",
    "aN":"2",
    "bN":"2",
    "cN":"2",
    "clust_rad":"5",
    "traAN":"2",
    "traBN":"2",
    "traCN":"2",
    "self_consistent":"0", # gets cast to int then bool
    "sc_temp":"sc_temp.template",
    "dev_tol":"0.001"}

    usr_inputs = rf.read_config(name)
    inputs.update(usr_inputs)

    return inputs
