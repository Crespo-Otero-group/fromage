import fromage.io.read_file as rf
from fromage.utils.dimer import Dimer

def test_init():
    big_mol = rf.mol_from_file("h2o_dimer.xyz")
    mol_a = big_mol.select(0)
    mol_b = big_mol.select(2)
    dim = Dimer([mol_a,mol_b])
    return

def test_dim_from_file():
    dim = rf.dim_from_file("h2o_dimer.xyz")
    print(dim)
