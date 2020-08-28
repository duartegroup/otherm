from .. import otherm
import pytest
import numpy as np
import os
here = os.path.dirname(os.path.abspath(__file__))

methane_path = os.path.join(here, 'data', 'methane.out')


def test_methane():

    methane = otherm.Molecule(methane_path)

    assert not methane.is_ts
    assert methane.n_atoms == 5
    assert len(methane.xyzs) == 5

    # Default rotational symmetry number is 1
    assert methane.sigma_r == 1

    assert methane.freqs is not None
    assert len(methane.freqs) == 3 * 5

    expected_freqs = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1260.03, 1261.06,
                      1263.09, 1483.73, 1483.97, 2967.29, 3109.79, 3111.75,
                      3112.88]

    # Ensure the frequencies have been extracted correctly and ordered
    # high to low
    for freq, expected in zip(methane.freqs[::-1], expected_freqs):
        assert np.abs(freq - expected) < 1E-6

    assert len(methane.real_vib_freqs()) == 3 * 5 - 6

    # And has the expected electronic energy
    expected_e = -40.416406478325
    expected_e_j = otherm.Constants.ha_to_j_mol * expected_e

    assert np.abs(methane.e - expected_e_j) < 1E-6

    # Carbon atom should be ~ at the origin once shifted to the COM
    for i, (atom_label, x, y, z) in enumerate(methane.xyzs):

        if atom_label != 'C':
            continue

        assert np.linalg.norm(methane.coords()[i]) < 1E-3

    # Total mass of methane is ~16 amu
    assert np.abs(methane.mass - otherm.Constants.amu_to_kg * 16.04) < 1E-3

    # Ensure the pcoords have the correct shape
    pcoords = methane.pcoords()
    assert pcoords.shape == (2, 5, 3)


def test_bad_molecule():

    with pytest.raises(Exception):
        _ = otherm.Molecule('not_a_file')

        _ = otherm.Molecule('README.md')
