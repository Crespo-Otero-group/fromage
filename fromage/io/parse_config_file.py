"""Reads the user's config inputs

This contains all of the default settings and may raise warnings and exceptions
if contradictory features are included.
"""
import numpy as np

from fromage.io import read_file as rf


def complete_config(name="config"):
    """Write default parameters for config and update with user inputs"""
    inputs = {
        "name": "fromage_calc",
        "target_shell": "",
        "cell_file": "cell.xyz",
        "high_pop_program": "gaussian",
        "high_pop_file": "gaussian_h.log",
        "high_pop_method": "ESP",
        "low_pop_program": "gaussian",
        "low_pop_file": "gaussian_l.log",
        "low_pop_method": "ESP",
        "bonding":"dis",
        "bond_thresh": "1.7",
        "atom_label": "1",
        "ewald": "",  # becomes bool
        "nchk": "1000",
        "nat": "500",
        "an": "2",
        "bn": "2",
        "cn": "2",
        "a_vec": "",
        "b_vec": "",
        "c_vec": "",
        "vectors_file": "",
        "clust_rad": "5",
        "clust_mode": "inc",
        "traan": "2",
        "trabn": "2",
        "tracn": "2",
        "self_consistent": "",  # becomes bool
        "sc_temp": "sc_temp.template",
        "dev_tol": "0.001",
        "damping": "0.0",
        "print_tweak": ""}

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

    inputs["bond_thresh"] = float(inputs["bond_thresh"])
    if type(inputs["atom_label"]) == str:
        inputs["atom_label"] = [int(inputs["atom_label"]) - 1]
    else:
        inputs["atom_label"] = [int(i) - 1 for i in inputs["atom_label"]]
    inputs["ewald"] = bool_cast(inputs["ewald"])
    inputs["nchk"] = int(inputs["nchk"])
    inputs["nat"] = int(inputs["nat"])
    inputs["an"] = int(inputs["an"])
    inputs["bn"] = int(inputs["bn"])
    inputs["cn"] = int(inputs["cn"])
    inputs["clust_rad"] = float(inputs["clust_rad"])
    inputs["traan"] = int(inputs["traan"])
    inputs["trabn"] = int(inputs["trabn"])
    inputs["tracn"] = int(inputs["tracn"])
    inputs["self_consistent"] = bool_cast(inputs["self_consistent"])
    inputs["dev_tol"] = float(inputs["dev_tol"])
    inputs["a_vec"] = np.array([float(i) for i in inputs["a_vec"]])
    inputs["b_vec"] = np.array([float(i) for i in inputs["b_vec"]])
    inputs["c_vec"] = np.array([float(i) for i in inputs["c_vec"]])
    inputs["damping"] = float(inputs["damping"])
    inputs["print_tweak"] = bool_cast(inputs["print_tweak"])
    # specified in config
    inputs["vectors"] = np.zeros((3, 3))
    if all(len(i) == 3 for i in [inputs["a_vec"], inputs["b_vec"], inputs["c_vec"]]):
        a_vec = inputs["a_vec"]
        b_vec = inputs["b_vec"]
        c_vec = inputs["c_vec"]
        inputs["vectors"][0] = a_vec
        inputs["vectors"][1] = b_vec
        inputs["vectors"][2] = c_vec
    else:  # from external file
        inputs["vectors"] = rf.read_vectors(inputs["vectors_file"])

    return inputs
