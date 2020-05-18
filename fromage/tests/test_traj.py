def test_read_traj_xyz(h2o_dimer_traj):
    assert len(h2o_dimer_traj) == 5

def test_read_traj_outcar(mbi_opt_traj):
    assert len(mbi_opt_traj) == 5

def test_write(h2o_dimer_traj):
    h2o_dimer_traj.write_xyz("boop.xyz")
    assert 1==1


