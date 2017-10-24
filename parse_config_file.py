#!/usr/bin/env python
"""Reads the user's config inputs

This contains all of the default settings and may raise warnings and exceptions
if contradictory features are included.
"""
import read_file as rf
import numpy as np

def complete_config(name="config"):
    """Write default parameters for config and update with user inputs"""
    inputs = {
        "cell_file": "cell.xyz",
        "high_pop_program": "gaussian",
        "high_gauss_file": "gaussian.log",
        "high_gauss_method": "ESP",
        "high_cp2k_file": "cp2k.out",
        "high_cp2k_method": "ESP",
        "low_pop_program": "gaussian",
        "low_gauss_file": "gaussian.log",
        "low_gauss_method": "ESP",
        "low_cp2k_file": "cp2k.out",
        "low_cp2k_method": "ESP",
        "max_bl": "1.7",
        "atom_label": "1",
        "ewald": "",  # gets cast to int then bool
        "nChk": "1000",
        "nAt": "500",
        "aN": "2",
        "bN": "2",
        "cN": "2",
        "clust_rad": "5",
        "traAN": "2",
        "traBN": "2",
        "traCN": "2",
        "self_consistent": "",  # gets cast to int then bool
        "sc_temp": "sc_temp.template",
        "dev_tol": "0.001"}

    usr_inputs = rf.read_config(name)
    inputs.update(usr_inputs)

    # Parse dictionary where necessary
    inputs["cell_file"]

    return inputs


def isfloat(in_str):
    """Check if a string can be cast to float"""
    try:
        float(in_str)
        return True
    except ValueError:
        return False

def bool_cast(in_str):
    """Casts string to bool avoiding some traps"""
    if in_str.lower().strip() in ("false","no","off","zero","none","","nan"):
        out_bool = False
    elif isfloat(in_str):
        out_bool = bool(float(in_str))
    else:
        out_bool = bool(in_str)
    return out_bool

def parse_inputs(name="config"):
    """Convert the string values of the diciontary to the appropriate types"""

    inputs=complete_config(name)
    inputs["max_bl"] = float(inputs["max_bl"])
    inputs["atom_label"] = int(inputs["atom_label"])
    inputs["ewald"] = bool_cast(inputs["max_bl"])
    inputs["nChk"] = int(inputs["nChk"])
    inputs["nAt"] = int(inputs["nAt"])
    inputs["aN"] = int(inputs["aN"])
    inputs["bN"] = int(inputs["bN"])
    inputs["cN"] = int(inputs["cN"])
    inputs["clust_rad"] = int(inputs["clust_rad"])
    inputs["traAN"] = int(inputs["traAN"])
    inputs["traBN"] = int(inputs["traBN"])
    inputs["traCN"] = int(inputs["traCN"])
    inputs["self_consistent"] = bool_cast(inputs["nChk"])
    inputs["dev_tol"] = int(inputs["dev_tol"])
    inputs["a_vec"]=np.array(inputs["a_vec"])
    inputs["b_vec"]=np.array(inputs["b_vec"])
    inputs["c_vec"]=np.array(inputs["c_vec"])

    return inputs
