def test_select_h2o_dimer_mol(h2o_dimer_mol):
    water = h2o_dimer_mol.select(3)
    assert len(water) == 3


def test_select_hc1_quad(hc1_quad):
    mol = hc1_quad.select([0])
    assert len(mol) == 37


def test_multiselect_hc1_quad(hc1_quad):
    mol = hc1_quad.select([0, 74])
    assert len(mol) == 74


def test_per_select_hc1_cell(hc1_cell):
    """Check that periodic select gets all atoms"""
    selected = hc1_cell.per_select(0)
    assert len(selected) == 37


def test_per_select_complete(hc1_cell, hc1_quad):
    """Check that periodic select completes the molecule"""
    selected = hc1_cell.per_select(0)
    new, old = hc1_cell.per_select(0, old_pos=True)
    new_sel = selected.select(0)
    assert len(selected) == len(new_sel)


def test_multi_per_select(hc1_cell, hc1_quad):
    selected = hc1_cell.per_select([0, 1])
    new_sel = hc1_quad.select(0)
    assert len(selected) == len(new_sel) * 2


def test_mol_segregation(hc1_quad):
    mols = hc1_quad.segregate()
    assert len(mols) == 4 and len(mols[1]) == 37
