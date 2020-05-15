def test_read_traj(h2o_dimer_traj):
    assert len(h2o_dimer_traj) == 5

def test_write(h2o_dimer_traj):
    h2o_dimer_traj.write_xyz("boop.xyz")
    assert 1==1


