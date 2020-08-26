from .. import otherm


def test_methane():

    methane = otherm.Molecule('methane.out')
    sigma_r = otherm.calc_symmetry_number(molecule=methane,
                                          max_n_fold_rot_searched=6,
                                          dist_tol=0.1)
    assert sigma_r == 12


def test_benzene():

    methane = otherm.Molecule('benzene.out')
    sigma_r = otherm.calc_symmetry_number(molecule=methane,
                                          max_n_fold_rot_searched=6,
                                          dist_tol=0.2)
    assert sigma_r == 12
