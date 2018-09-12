import pytest
from pytest import approx
import cryspy.io.read_file as rf
import numpy as np

@pytest.fixture
def benz_cell():
    """Return a Mol object of the benzene cell"""
    out_cell = rf.mol_from_file("benzene_cell.xyz")
    out_cell.vectors = rf.read_vectors("benzene_vectors")
    return out_cell

def test_import(benz_cell):
    print()
    print(benz_cell)
    print("hello")
    return
