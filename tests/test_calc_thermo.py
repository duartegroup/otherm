from .. import otherm
import numpy as np


def test_methane():

    methane = otherm.Molecule('methane.out')

    # ORCA 4.0 defaulted to sigma_r = 4 and 1 atmosphere standard state
    methane.calculate_thermochemistry(ss='1atm',
                                      symm_n=3)

    expected_g = -40.39166850 * otherm.Constants.ha_to_j_mol

    # Difference should be less than 1 J mol-1
    assert np.abs(methane.g - expected_g) < 1


def test_methane_gaussian():

    methane = otherm.Molecule('methane.out')

    # Populate parameters from a Gaussian09 calculation
    methane.freqs = (6 * [0.0] +
                     [1308.0624, 1308.0891, 1308.1077, 1530.3570, 1530.3806,
                      3046.1934, 3198.8181, 3198.8521, 3198.9170])[::-1]

    methane.xyzs = [['C', -0.000002,  0.000041, 0.000008],
                    ['H', -0.637549, -0.828986, -0.334635],
                    ['H', -0.444245, 0.953621, -0.314719],
                    ['H', 0.083441, -0.021253,  1.094689],
                    ['H', 0.998355, -0.103324, -0.445342]]
    methane.shift_to_com()

    expected_i_mat = np.array([[-0.34825, 0.18296, 0.91937],
                               [0.67773, 0.72672, 0.11209],
                               [0.64762, -0.66212, 0.37707]])

    i_mat = otherm.calc_moments_of_inertia(methane.xyzs)
    assert np.sum(i_mat - expected_i_mat) < 1E-6

    methane.e = -40.4294662 * otherm.Constants.ha_to_j_mol
    expected_g = -40.404434 * otherm.Constants.ha_to_j_mol

    # Thermochemistry not calculated with symmetry
    methane.calculate_thermochemistry(ss='1atm',
                                      method='igm',
                                      symm_n=1)

    # Slightly poorer agreement with Gaussian but still within 0.5 kJ mol-1
    assert np.abs(methane.g - expected_g) < 4


def test_calc_ss():

    methane_1atm = otherm.Molecule('methane.out')
    methane_1atm.calculate_thermochemistry(ss='1atm')

    methane_1m = otherm.Molecule('methane.out')
    methane_1m.calculate_thermochemistry(ss='1M')

    delta = methane_1m.g - methane_1atm.g
    # Expected additive amount it 1.9 kcal mol-1

    assert 1.8 < otherm.Constants.j_to_kcal * delta < 2
