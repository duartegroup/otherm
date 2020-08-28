from .. import otherm
import numpy as np
import os
here = os.path.dirname(os.path.abspath(__file__))

methane_path = os.path.join(here, 'data', 'methane.out')
benzene_path = os.path.join(here, 'data', 'benzene.out')


def path(out_name):
    return os.path.join(here, 'data', out_name)


def test_same_under():

    methane = otherm.Molecule(methane_path)
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


def test_highest_cn():

    methane = otherm.Molecule(methane_path)

    cn_axes = otherm.cn_and_axes(methane, max_n=6, dist_tol=0.25)
    # 3 C2 axes
    assert len(cn_axes[2]) == 3

    # Should have 4 C3 axes along each CH bond
    assert len(cn_axes[3]) == 4

    benzene = otherm.Molecule(benzene_path)
    n_c6 = len(otherm.cn_and_axes(benzene, max_n=8, dist_tol=0.25)[6])
    # benzene should have a single C6 axis
    assert n_c6 == 1

    # Expected number of C2, C3, ... C8 axes in benzene. Includes e.g C6^3 = C2
    expected_n_axes = {2: 7, 3: 1, 4: 0, 5: 0, 6: 1, 7: 0, 8: 0}

    for n, axes in otherm.cn_and_axes(benzene, max_n=8, dist_tol=0.25).items():
        assert expected_n_axes[n] == len(axes)


def test_symmetry_number():

    methane = otherm.Molecule(path('methane.out'))

    assert not methane.is_linear()
    assert otherm.calc_symmetry_number(methane) == 12

    benzene = otherm.Molecule(path('benzene.out'))
    assert otherm.calc_symmetry_number(benzene) == 12

    bh3 = otherm.Molecule(path('BH3.out'))
    assert otherm.calc_symmetry_number(bh3) == 6       # D3h

    co = otherm.Molecule(path('CO.out'))
    assert otherm.calc_symmetry_number(co) == 1        # Cinf_v

    co2 = otherm.Molecule(path('CO2.out'))
    assert otherm.calc_symmetry_number(co2) == 2       # Dinf_h

    cp = otherm.Molecule(path('cp.out'))
    assert otherm.calc_symmetry_number(cp) == 10       # D5h

    ethane = otherm.Molecule(path('ethane.out'))
    assert otherm.calc_symmetry_number(ethane) == 6    # D3d

    ethene = otherm.Molecule(path('ethene.out'))
    assert otherm.calc_symmetry_number(ethene) == 4    # D2h

    h2o = otherm.Molecule(path('H2O.out'))
    assert otherm.calc_symmetry_number(h2o) == 2       # C2v

    nh3 = otherm.Molecule(path('nh3.out'))
    assert otherm.calc_symmetry_number(nh3) == 3       # C3v
