from pytest import approx
import numpy as np
from fromage.utils.atom import Atom
from fromage.utils.mol import Mol


def test_at_list_type(at_list):
    """The atom list is made of Atom objects"""
    assert all(isinstance(i, Atom) for i in at_list)


def test_init(at_list):
    """The Mol object is initialised correctly"""
    mo = Mol(at_list)
    assert isinstance(mo, Mol)


def test_duplicity(h2o_dup):
    """Check for repeated atoms"""
    h2o_dup.remove_duplicates()
    assert len(h2o_dup) == 3


def test_len(h2o_dimer_mol):
    """The len method is implemented"""
    assert len(h2o_dimer_mol) == 6


def test_len_empty(empty_mol):
    assert len(empty_mol) == 0


def test_for_loop(h2o_dimer_mol):
    counter = 0
    for at in h2o_dimer_mol:
        counter += 1
    assert counter == 6


def test_add(h2o_dimer_mol, newat):
    """Addition is implemented"""
    new_mol = h2o_dimer_mol + Mol(newat)
    assert len(new_mol) == 7


def test_append(h2o_dimer_mol, newat):
    """Appending is implemented"""
    h2o_dimer_mol.append(newat)
    assert len(h2o_dimer_mol) == 7


def test_centroid(c_o):
    cen = c_o.centroid()
    assert cen[0] == approx(0.5)


def test_center(c_o):
    c_o.center_mol()
    assert c_o[0].x == approx(-0.5)


def test_translate(c_o):
    vec = np.array([0.25, 1.0, 0.0])
    c_o.translate(vec)
    assert c_o[1].x == approx(1.25)


def test_set_bonding_str(h2o_dimer_mol):
    h2o_dimer_mol.set_bonding_str("dis")
    assert h2o_dimer_mol.bonding == "dis"
    assert h2o_dimer_mol.thresh == 1.8
    h2o_dimer_mol.set_bonding_str("1.5")
    assert h2o_dimer_mol.bonding == "dis"
    assert h2o_dimer_mol.thresh == 1.5
    h2o_dimer_mol.set_bonding_str("vdw")
    assert h2o_dimer_mol.bonding == "vdw"
    assert h2o_dimer_mol.thresh == -0.3
    h2o_dimer_mol.set_bonding_str("cov-0.1")
    assert h2o_dimer_mol.bonding == "cov"
    assert h2o_dimer_mol.thresh == -0.1
    h2o_dimer_mol.set_bonding_str("-0.1cov")
    assert h2o_dimer_mol.bonding == "cov"
    assert h2o_dimer_mol.thresh == -0.1

def test_same_atoms_as():
    mol_a_lis = [Atom("C",1.,0.,0.),
                 Atom("H",1.,1.,1.)]
    mol_b_lis = [Atom("H",1.,1.,1.),
                 Atom("C",1.,0.,0.)]
    mol_c_lis = [Atom("C",1.,0.,0.),
                 Atom("C",1.,1.,1.)]
    mol_a = Mol(mol_a_lis)
    mol_b = Mol(mol_b_lis)
    mol_c = Mol(mol_c_lis)

    assert mol_a.same_atoms_as(mol_b) == True
    assert mol_a.same_atoms_as(mol_c) == False
