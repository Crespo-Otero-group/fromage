import fromage.io.read_file as rf
from fromage.utils.dimer import Dimer
from fromage.tests.conftest import _in_data

def test_init(h2o_dimer):
    mol_a = h2o_dimer.select(0)
    mol_b = h2o_dimer.select(2)
    dim = Dimer([mol_a, mol_b])
    return


def test_dim_from_file():
    dim = rf.dim_from_file(_in_data("h2o_dimer.xyz"))
    print(dim)
