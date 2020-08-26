from .. import otherm
import numpy as np


def test_same_under1():

    methane = otherm.Molecule('methane.out')
    coords = methane.coords()

    # A C-H vector should be the final atom to another
    assert methane.xyzs[4][0] == 'C'
    assert methane.xyzs[3][0] == 'H'

    axis = coords[4] - coords[3]
    axis /= np.linalg.norm(axis)

    axes = otherm.get_possible_axes(methane.coords())
    assert any(np.linalg.norm(ax - axis) < 1E-3 for ax in axes)

    assert otherm.is_same_under_n_fold(pcoords=methane.pcoords(),
                                       axis=axis,
                                       n=3)

    # Should not be identical under a 4-fold rotation along a C-H bond
    assert not otherm.is_same_under_n_fold(pcoords=methane.pcoords(),
                                           axis=axis,
                                           n=4)

    # Should havea C2 along an axis pointing between two CH bonds
    ch1 = coords[3] - coords[4]
    ch2 = coords[2] - coords[4]

    axis = np.average(np.array([ch1, ch2]), axis=0)
    axis /= np.linalg.norm(axis)

    assert otherm.is_same_under_n_fold(pcoords=methane.pcoords(),
                                       axis=axis,
                                       n=2)


def test_is_same_under2():

    benzene = otherm.Molecule('benzene.out')
    # TODO: finish this test


def test_highest_cn():

    methane = otherm.Molecule('methane.out')

    cn_axes = otherm.cn_and_axes(methane, max_n=6, dist_tol=0.1)
    # 3 C2 axes
    assert len(cn_axes[2]) == 3

    # Should have 4 C3 axes along each CH bond
    assert len(cn_axes[3]) == 4

    benzene = otherm.Molecule('benzene.out')
    n_c6 = len(otherm.cn_and_axes(benzene, max_n=8, dist_tol=0.2)[6])
    # benzene should have a single C6 axis
    assert n_c6 == 1

    # Expected number of C2, C3, ... C8 axes in benzene. Includes e.g C6^3 = C2
    expected_n_axes = {2: 7, 3: 1, 4: 0, 5: 0, 6: 1, 7: 0, 8: 0}

    for n, axes in otherm.cn_and_axes(benzene, max_n=8, dist_tol=0.25).items():
        assert expected_n_axes[n] == len(axes)


def test_symmetry_number():

    methane = otherm.Molecule('methane.out')

    assert not otherm.is_linear(coords=methane.coords())
    assert otherm.calc_symmetry_number(methane, dist_tol=0.1) == 12

    benzene = otherm.Molecule('benzene.out')
    assert otherm.calc_symmetry_number(benzene, dist_tol=0.25) == 12
