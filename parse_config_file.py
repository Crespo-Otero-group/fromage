"""Reads the user's config inputs

This contains all of the default settings and may raise warnings and exceptions
if contradictory features are included.
"""
import read_file as rf
import numpy as np


def complete_config(name="config"):
    """Write default parameters for config and update with user inputs"""
    inputs = {
        "name":"cryspy_calc",
        "cell_file": "cell.xyz",
        "high_pop_program": "gaussian",
        "high_pop_file": "cp2k.out",
        "high_pop_method": "ESP",
        "low_pop_program": "gaussian",
        "low_pop_file": "gaussian.log",
        "low_pop_method": "ESP",
        "max_bl": "1.7",
        "atom_label": "1",
        "ewald": "",  # gets cast to int then bool
        "nchk": "1000",
        "nat": "500",
        "an": "2",
        "bn": "2",
        "cn": "2",
        "clust_rad": "5",
        "traan": "2",
        "trabn": "2",
        "tracn": "2",
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
    if in_str.lower().strip() in ("false", "no", "off", "zero", "none", "", "nan"):
        out_bool = False
    elif isfloat(in_str):
        out_bool = bool(float(in_str))
    else:
        out_bool = bool(in_str)
    return out_bool


def parse_inputs(name="config"):
    """Convert the string values of the diciontary to the appropriate types"""

    inputs = complete_config(name)

    inputs["max_bl"] = float(inputs["max_bl"])
    if type(inputs["atom_label"]):
        inputs["atom_label"] = [int(inputs["atom_label"]) - 1]
    else:
        label_atom = [int(i) - 1 for i in inputs["atom_label"]]
    inputs["ewald"] = bool_cast(inputs["ewald"])
    inputs["nchk"] = int(inputs["nchk"])
    inputs["nat"] = int(inputs["nat"])
    inputs["an"] = int(inputs["an"])
    inputs["bn"] = int(inputs["bn"])
    inputs["cn"] = int(inputs["cn"])
    inputs["clust_rad"] = int(inputs["clust_rad"])
    inputs["traan"] = int(inputs["traan"])
    inputs["trabn"] = int(inputs["trabn"])
    inputs["tracn"] = int(inputs["tracn"])
    inputs["self_consistent"] = bool_cast(inputs["self_consistent"])
    inputs["dev_tol"] = float(inputs["dev_tol"])
    inputs["a_vec"] = np.array(inputs["a_vec"])
    inputs["b_vec"] = np.array(inputs["b_vec"])
    inputs["c_vec"] = np.array(inputs["c_vec"])

    return inputs
