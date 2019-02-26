import pytest
import fromage.io.read_file as rf
from pytest import approx
from fromage.utils.dimer import Dimer

def test_init():
    big_mol = rf.mol_from_file("h2o_dimer.xyz")
    mol_a = big_mol.select(0)
    mol_b = big_mol.select(2)
    print(mol_a)
    print(mol_b)

    dim = Dimer([mol_a,mol_b])
    print(dim)
    return
