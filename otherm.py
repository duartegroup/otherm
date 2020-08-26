#!/usr/bin/env python3
import argparse
import numpy as np
from copy import deepcopy


class Constants:
    k_b = 1.38064852E-23                # J K-1
    h = 6.62607004E-34                  # J s
    n_a = 6.022140857E23                # molecules mol-1
    ha_to_j = 4.359744650E-18           # J Ha-1
    ha_to_j_mol = ha_to_j * n_a         # J mol-1 Ha-1
    atm_to_pa = 101325                  # Pa
    dm_to_m = 0.1                       # m
    amu_to_kg = 1.660539040E-27         # Kg
    r = k_b * n_a                       # J K-1 mol-1
    c = 299792458                       # m s-1
    c_in_cm = c * 100                   # cm s-1
    ang_to_m = 1E-10                    # m

    #  Atomic weights in amu from:
    #  IUPAC-CIAWW's Atomic weights of the elements: Review 2000
    atomic_masses = {"H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182,
                     "B": 10.811, "C": 12.0107, "N": 14.0067, "O": 15.9994,
                     "F": 18.9984032, "Ne": 2.01797, "Na": 22.989770,
                     "Mg": 24.3050, "Al": 26.981538, "Si": 28.0855,
                     "P": 30.973761, "S": 32.065, "Cl": 35.453, "Ar": 39.948,
                     "K": 39.0983, "Ca": 40.078, "Sc": 44.955910, "Ti": 47.867,
                     "V": 50.9415, "Cr": 51.9961, "Mn": 54.938049,
                     "Fe": 55.845, "Co": 58.933200, "Ni": 58.6934,
                     "Cu": 63.546, "Zn": 65.409, "Ga": 69.723, "Ge": 72.64,
                     "As": 74.92160, "Se": 78.96, "Br": 79.904, "Kr": 83.798,
                     "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585, "Zr": 91.224,
                     "Nb": 92.90638, "Mo": 95.94, "Ru": 101.07,
                     "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682,
                     "Cd": 112.411, "In": 114.818, "Sn": 118.710,
                     "Sb": 121.760, "Te": 127.60, "I": 126.90447,
                     "Xe": 131.293, "Cs": 132.90545, "Ba": 137.327,
                     "La": 138.9055, "Ce": 140.116, "Pr": 140.90765,
                     "Nd": 144.24, "Sm": 150.36, "Eu": 151.964, "Gd": 157.25,
                     "Tb": 158.92534, "Dy": 162.500, "Ho": 164.93032,
                     "Er": 167.259, "Tm": 168.93421, "Yb": 173.04,
                     "Lu": 174.967, "Hf": 178.49, "Ta": 180.9479,
                     "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
                     "Pt": 195.078, "Au": 196.96655, "Hg": 200.59,
                     "Tl": 204.3833, "Pb": 207.2, "Bi": 208.98038,
                     "Th": 232.0381, "Pa": 231.03588, "U": 238.02891}


def get_args():
    """Get command line arguments with argparse"""

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", action='store',
                        help='.out file with freq calculation performed')

    parser.add_argument('-t', '--temp', type=float, default=298,
                        help="Temperature (K) at which to calculate G, H and "
                             "S. Default: %(default)s")

    parser.add_argument('-ss', '--standard_state', type=str, default='1M',
                        help="Standard state. 1atm for gas phase and 1M for "
                             "solution phase. Default: %(default)s")

    parser.add_argument('-m', '--method', default='grimme', nargs='?',
                        choices=['igm', 'truhlar', 'grimme'],
                        help='Method by which to calculate G, H and S. '
                             'Default: %(default)s')

    parser.add_argument('-s', '--shift', type=float, default=100,
                        help="Cut-off (in cm-1) to use in Truhlar's method. "
                             "Frequencies below this will be shifted to this "
                             "value. Default: %(default)s")

    parser.add_argument('-w', '--w0', type=float, default=100,
                        help="Value of w0 to use in Grimme's interpolation "
                             "method Chem. Eur. J. 2012, 18, 9955 eqn. 8. "
                             "Default: %(default)s")

    parser.add_argument('-a', '--alpha', type=float, default=4,
                        help="Value of alpha to use in Grimme's interpolation "
                             "method Chem. Eur. J. 2012, 18, 9955 eqn. 8. "
                             "Default: %(default)s")

    parser.add_argument('-cs', '--calc_sym', action='store_true', default=False,
                        help="Force calculation of symmetry number "
                             "(N^3 algorithm) used for n_atoms < 50."
                             " Default: %(default)s")

    parser.add_argument('-sn', '--symn', type=int, default=None,
                        help="Override the symmetry number calculation. "
                             "Default: %(default)s")

    parser.add_argument('-r', '--real_freqs', action='store_true', default=False,
                        help='Convert any imaginary frequencies to their real '
                             'counterparts')

    parser.add_argument('-ts', '--transition_state', action='store_true',
                        default=False,
                        help='This species is a transition state so the lowest'
                             'imaginary should be ignored in calculating the'
                             'thermochemical contributions')

    return parser.parse_args()


def print_output(molecule):

    print("----------------------------------------------------------------------------------\n"
          "|                                                                                 |\n"
          "|          /$$$$$$  /$$$$$$$$ /$$   /$$ /$$$$$$$$ /$$$$$$$  /$$      /$$          |\n"
          "|         /$$__  $$|__  $$__/| $$  | $$| $$_____/| $$__  $$| $$$    /$$$          |\n"
          "|        | $$  \ $$   | $$   | $$  | $$| $$      | $$  \ $$| $$$$  /$$$$          |\n"
          "|        | $$  | $$   | $$   | $$$$$$$$| $$$$$   | $$$$$$$/| $$ $$/$$ $$          |\n"
          "|        | $$  | $$   | $$   | $$__  $$| $$__/   | $$__  $$| $$  $$$| $$          |\n"
          "|        | $$  | $$   | $$   | $$  | $$| $$      | $$  \ $$| $$\  $ | $$          |\n"
          "|        |  $$$$$$/   | $$   | $$  | $$| $$$$$$$$| $$  | $$| $$ \/  | $$          |\n"
          "|         \______/    |__/   |__/  |__/|________/|__/  |__/|__/     |__/          |\n"
          "|                                                                                 |\n"
          "-----------------------------------------------------------------------------------\n\n")

    print("{:<50s}{:>33s}".format('Filename', args.filename))
    print("{:<50s}{:>33.1f}".format('Temperature (K)', args.temp))
    print("{:<50s}{:>33s}".format('Standard state is', args.standard_state))
    if args.real_freqs:
        print("{:<50s}{:>33s}".format('Treat imaginary (negative) frequencies '
                                      'as real', 'True'))

    print("{:<50s}{:>33s}".format('Calculating using the method of',
                                  args.method))
    if args.method.lower() == 'grimme':
        print("{:<50s}{:>33s}".format('', 'Chem. Eur. J. 2012, 18, 9955'))

    if args.method.lower() == 'truhlar':
        print("{:<50s}{:>33s}".format('', 'J. Phys. Chem. B, 2011, 115, 14556'))

    print()
    print("{:<50s}{:>33s}".format('Symmetry number (σ)', str(molecule.sigma_r)))
    print("{:<50s}{:>33.2f}".format('Molecular weight (amu)', molecule.mass / Constants.amu_to_kg))
    print()
    if any(freq < 0 for freq in molecule.freqs):
        print('---------------------------------------WARNING--------------------------------------')
        print('                   Found imaginary frequencies', )
        print('------------------------------------------------------------------------------------\n')
    print("{:<50s}{:>33.2f}".format('Total entropy (J K-1 mol-1)', molecule.s))
    print("{:<50s}{:>33.2f}".format('Total enthalpy (J mol-1)', molecule.h))
    print("{:<50s}{:>33.2f}".format('Total free energy (J mol-1)', molecule.g))
    print()

    e_elec_ha = molecule.e / (Constants.ha_to_j * Constants.n_a)
    h_ha = molecule.h / (Constants.ha_to_j * Constants.n_a)
    g_ha = molecule.g / (Constants.ha_to_j * Constants.n_a)

    print("{:<50s}{:>33}".format('For convenience E, H, G in Hartrees', ''))
    print(e_elec_ha, h_ha, g_ha, sep=',')
    print("----------------------------------------------------------------------------------")

    return None


def extract_frequencies(filename):
    """
    Extract the frequencies from an ORCA output file. Iterate through the
    reversed file to find the frequencies (in cm-1) so if multiple Hessians
    have been calculated, the final one is found.

    :param filename: Name of the ORCA output file
    :return: (list(float)) List of frequencies (high to low, in cm-1)
    """
    orca_out_file_lines = open(filename, 'r').readlines()

    freq_block = False
    freq_list = []

    for line in reversed(orca_out_file_lines):

        if 'NORMAL MODES' in line:
            freq_block = True
        if 'VIBRATIONAL FREQUENCIES' in line:
            break
        if 'cm**-1' in line and freq_block:
            try:
                # Line is in the form "   0:         0.00 cm**-1"
                freq_list.append(float(line.split()[1]))

            except TypeError:
                raise Exception("Couldn't extract frequencies")

    if len(freq_list) == 0:
        raise Exception('Frequencies not found')

    return freq_list


def extract_final_electronic_energy(filename):
    """
    Get the final electronic energy from an ORCA output file

    :param filename: (str)
    :return: (float)
    """

    for line in reversed(open(filename, 'r').readlines()):
        if 'FINAL SINGLE POINT ENERGY' in line:
            return Constants.ha_to_j * Constants.n_a * float(line.split()[4])

    raise Exception('Final electronic energy not found')


def extract_xyzs(filename):
    """
    Extract the xyzs from a ORCA output file in the format [[C, 0.0000, 0.0000, 0.0000], ....]
    :param filename: Name of the ORCA output file
    :return: List of xyzs
    """

    xyzs = []
    orca_out_file_lines = []

    try:
        orca_out_file_lines = [line for line in open(filename, 'r')]
    except IOError as err:
        exit("I/O error({0}): {1}".format(err.errno, err.strerror))

    cart_coords_block = False

    for line in reversed(orca_out_file_lines):

        xyz_block = True if len(line.split()) == 4 else False

        if cart_coords_block and xyz_block:
            atom, x, y, z = line.split()[-4:]
            xyzs.append([atom, float(x), float(y), float(z)])

        if 'CARTESIAN COORDINATES (A.U.)' in line:
            cart_coords_block = True

        if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
            break

    return xyzs


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def is_linear(coords, cos_angle_tol=0.05):
    """
    Determine if a molecule is linear

    :param coords:
    :param cos_angle_tol:
    :return:
    """

    if len(coords) == 2:
        return True

    if len(coords) > 2:
        vec_atom0_atom1 = calc_normalised_vector(coords[0], coords[1])
        is_atom_colinear = [False for _ in range(2, len(coords))]
        for i in range(2, len(coords)):
            vec_atom0_atomi = calc_normalised_vector(coords[0], coords[i])
            if 1.0 - cos_angle_tol < np.abs(np.dot(vec_atom0_atom1, vec_atom0_atomi)) < 1.0:
                is_atom_colinear[i - 2] = True

        if all(is_atom_colinear):
            return True

    return False


def calc_normalised_vector(coord1, coord2):
    vec = coord2 - coord1
    return vec / np.linalg.norm(vec)


def is_same_structure(old_coords, new_coords, xyzs, dist_tol):
    """
    Determine whether a the new structure maps onto the old, i.e is the same structure under some symmetry operation

    :param old_coords:
    :param new_coords:
    :param xyzs:
    :param dist_tol:
    :return:
    """

    n_atoms = len(xyzs)
    close_pairs = []

    for j in range(n_atoms):
        for k in range(n_atoms):
            if xyzs[j][0] == xyzs[k][0]:
                if np.linalg.norm(new_coords[j] - old_coords[k]) < dist_tol:
                    close_pairs.append([j, k])

        if len(close_pairs) != j + 1:                       # If atom j doesn't map then we can break immediately
            return False

    if len(close_pairs) == n_atoms:
        return True

    return False


def strip_identical_and_inv_axes(axes, sim_axis_tol):
    """
    For a list of axes remove those which are similar to within some distance tolerance, or are inverses to within
    that tolerance

    :param axes: list of axes
    :param sim_axis_tol: distance tolerance in Å
    :return:
    """

    unique_possible_axes = []

    for i in range(len(axes)):
        unique = True
        for unique_axis in unique_possible_axes:
            if np.linalg.norm(axes[i] - unique_axis) < sim_axis_tol:
                unique = False
            if np.linalg.norm(-axes[i] - unique_axis) < sim_axis_tol:
                unique = False
        if unique:
            unique_possible_axes.append(axes[i])

    return unique_possible_axes


def get_possible_axes(coords, max_triple_dist=2.0, sim_axis_tol=0.05):
    """
    Possible rotation axes in a molecule. Currently limited to average vectors and cross products i.e.

      Y          Y --->
     / \        / \
    X   Y      X   Z

      |
      |
      ,

    :param coords:
    :return:
    """

    possible_axes = []
    n_atoms = len(coords)

    for i in range(n_atoms):
        for j in range(n_atoms):

            if i > j:                   # For the unique pairs add the i–j vector
                possible_axes.append(calc_normalised_vector(coords[i], coords[j]))

            for k in range(n_atoms):
                if i != j and i != k and j != k:
                    vec1 = coords[j] - coords[i]
                    vec2 = coords[k] - coords[i]
                    if np.linalg.norm(vec1) < max_triple_dist and np.linalg.norm(vec2) < max_triple_dist:
                            avg_vec = (vec1 + vec2) / 2.0
                            possible_axes.append(avg_vec/np.linalg.norm(avg_vec))
                            perp_vec = np.cross(vec1, vec2)
                            possible_axes.append(perp_vec/np.linalg.norm(perp_vec))

    unique_possible_axes = strip_identical_and_inv_axes(possible_axes, sim_axis_tol)

    return unique_possible_axes


def apply_n_fold_rotation(coords, axis, n):
    """
    This function will SHIFT the origin of the coordinates to (0.0, 0.0, 0.0)

    :param coords:
    :param origin:
    :param axis:
    :param n:
    :return:
    """

    tmp_coords = np.array(coords, copy=True)
    rot_mat = rotation_matrix(axis, theta=(2.0 * np.pi / n))
    return np.array([np.matmul(rot_mat, coord) for coord in tmp_coords])


def find_highest_cn(coords, xyzs, max_n, dist_tol):
    """
    Find the highest symmetry rotation axis

    :param coords:
    :param xyzs:
    :param max_n:
    :param dist_tol:
    :return:
    """
    axes = get_possible_axes(coords)

    n_axes = [[] for _ in range(max_n+1)]

    for axis in axes:

        # Minimum n-fold rotation is 2
        for n in range(2, max_n+1):
            tmp_coords = apply_n_fold_rotation(coords, axis, n)

            if is_same_structure(tmp_coords, coords, xyzs, dist_tol):
                n_axes[n].append(axis)

    max_n, axes_max_n = 0, []
    for n, axes in enumerate(n_axes):
        if len(axes) > 0:
            axes_max_n, max_n = axes, n

    return max_n, axes_max_n


def calc_symmetry_number(molecule, max_n_fold_rot_searched=6, dist_tol=0.2):
    """
    Calculate the symmetry number of a molecule.
    Based on Theor Chem Account (2007) 118:813–826

    :param molecule:
    :param max_n_fold_rot_searched:
    :param dist_tol:
    :return:
    """
    symmetry_number = 1
    coords = molecule.coords()

    # Ensure the origin is at the center of mass
    if np.max(molecule.com) > 0.1:
        molecule.shift_to_com()

    max_n, axes = find_highest_cn(coords,
                                  molecule.xyzs,
                                  max_n_fold_rot_searched, dist_tol)
    if max_n == 1:
        # C1 symmetry → σ=1
        return 1

    possible_axes = get_possible_axes(coords)
    curr_symmetry_number = 1                          # Already has E symmetry

    # For every possible axis apply C2, C3...C_n_max rotations
    for axis in possible_axes:
        for n in range(2, max_n_fold_rot_searched+1):

            tmp_coords = apply_n_fold_rotation(coords, axis, n)

            if is_same_structure(tmp_coords, coords, molecule.xyzs, dist_tol):
                if n == 2:
                    curr_symmetry_number += 1
                if n > 2:
                    curr_symmetry_number += 2

    if curr_symmetry_number > symmetry_number:
        symmetry_number = curr_symmetry_number

    if is_linear(coords):          # If the molecule is linear the then the symmetry_number will be wrong... fix it here
        if symmetry_number > 6:    # There are perpendicular C2s the point group is D∞h
            return 2
        else:                      # If not then C∞v and the symmetry number is 1
            return 1

    return symmetry_number


def calc_moments_of_inertia(xyz_list):
    """
    From a list of xyzs compute the matrix of moments of inertia
    :param xyz_list: List of xyzs in the format [[C, 0.0000, 0.0000, 0.0000], ....]
    :return: The matrix
    """

    i_matrix = np.zeros([3, 3])

    for xyz_line in xyz_list:
        atom_label, x_ang, y_ang, z_ang = xyz_line
        x, y, z = Constants.ang_to_m * x_ang, Constants.ang_to_m * y_ang, Constants.ang_to_m * z_ang
        atom_mass_kg = Constants.atomic_masses[atom_label] * Constants.amu_to_kg

        i_matrix[0, 0] += atom_mass_kg * (y**2 + z**2)
        i_matrix[0, 1] += -atom_mass_kg * (x * y)
        i_matrix[0, 2] += -atom_mass_kg * (x * z)

        i_matrix[1, 0] += -atom_mass_kg * (y * x)
        i_matrix[1, 1] += atom_mass_kg * (x**2 + z**2)
        i_matrix[1, 2] += -atom_mass_kg * (y * z)

        i_matrix[2, 0] += -atom_mass_kg * (z * x)
        i_matrix[2, 1] += -atom_mass_kg * (z * y)
        i_matrix[2, 2] += atom_mass_kg * (x**2 + y**2)

    return i_matrix


def calc_q_trans_igm(molecule, ss, temp):
    """
    Calculate the translational partition function using the PIB model,
    coupled with an effective volume

    :param molecule: (otherm.Molecule)
    :param ss: (str) Standard state to use {1M, 1atm}
    :param temp: (float) Temperature in K
    :return: (float) Translational partition function q_trns
    """

    if ss.lower() == '1atm':
        effective_volume = Constants.k_b * temp / Constants.atm_to_pa

    elif ss.lower() == '1m':
        effective_volume = 1.0 / (Constants.n_a * (1.0 / Constants.dm_to_m)**3)

    else:
        raise NotImplementedError

    molecule.q_trans = ((2.0 * np.pi * molecule.mass * Constants.k_b * temp / Constants.h**2)**1.5 *
                        effective_volume)

    return molecule.q_trans


def calc_q_rot_igm(molecule, temp):
    """
    Calculate the rotational partition function using the IGM method. Uses the
    rotational symmetry number, default = 1

    :param molecule: (otherm.Molecule)
    :param temp: (float) Temperature in K
    :return: (float) Rotational partition function q_rot
    """

    i_mat = calc_moments_of_inertia(molecule.xyzs)

    omega_x = Constants.h**2 / (8.0 * np.pi**2 * Constants.k_b * i_mat[0, 0])
    omega_y = Constants.h**2 / (8.0 * np.pi**2 * Constants.k_b * i_mat[1, 1])
    omega_z = Constants.h**2 / (8.0 * np.pi**2 * Constants.k_b * i_mat[2, 2])

    if molecule.n_atoms == 1:
        molecule.q_rot = 1
        return molecule.q_rot

    else:
        molecule.q_rot = (np.sqrt(np.pi)/molecule.sigma_r) * (temp**1.5 / np.sqrt(omega_x * omega_y * omega_z))

    return molecule.q_rot


def calc_q_vib_igm(molecule, temp):
    """
    Calculate the vibrational partition function using the IGM method.
    Uses the rotational symmetry number, default = 1

    :param molecule: (otherm.Molecule)
    :param temp: (float) Temperature in K
    :return: (float) Vibrational partition function q_rot
    """

    if molecule.n_atoms == 1:
        molecule.q_vib = 1
        return molecule.q_vib

    for freq in molecule.real_vib_freqs():
        x = freq * Constants.c_in_cm * Constants.h / Constants.k_b
        molecule.q_vib *= np.exp(-x / (2.0 * temp)) / (1.0 - np.exp(-x / temp))

    return molecule.q_vib


def calc_s_trans_pib(molecule, ss, temp):
    """Calculate the translational entropy using a particle in a box model

    :param molecule: (otherm.Molecule)
    :param ss: (str) Standard state to use for calculating the effective box
               size in the q_trans calculation
    :param temp: (float) Temperature in K
    :return: (float) S_trans in J K-1 mol-1
    """

    q_trans = calc_q_trans_igm(molecule, ss=ss, temp=temp)
    return Constants.r * (np.log(q_trans) + 1.0 + 1.5)


def calc_s_rot_rr(molecule, temp):
    """
    Calculate the rigid rotor (RR) entropy

    :param molecule: (otherm.Molecule)
    :return: (float) S_rot in J K-1 mol-1
    """

    if molecule.n_atoms == 1:
        return 0

    q_rot = calc_q_rot_igm(molecule, temp=temp)
    return Constants.r * (np.log(q_rot) + 1.5)


def calc_igm_s_vib(molecule, temp):
    """
    Calculate the entropy of a molecule according to the Ideal Gas Model (IGM)
    RRHO method

    :param molecule: (otherm.Molecule)
    :param temp: (float) Temperature in K
    :return: (float) S_vib in J K-1 mol-1
    """
    s = 0.0

    for freq in molecule.real_vib_freqs():
        x = freq * Constants.c_in_cm * Constants.h / (Constants.k_b * temp)
        s += Constants.r * ((x / (np.exp(x) - 1.0)) - np.log(1.0 - np.exp(-x)))

    return s


def calc_truhlar_s_vib(molecule, temp, shift_freq):
    """
    Calculate the entropy of a molecule according to the Truhlar's method of
    shifting low frequency modes

    :param molecule: (otherm.Molecule)
    :param temp: (float) Temperature in K
    :param shift_freq: (float) Shift all frequencies to this value
    :return: (float) S_vib in J K-1 mol-1
    """
    s = 0

    for freq in molecule.real_vib_freqs():

        # Threshold lower bound of the frequency, includes *all* imaginary
        freq = max(freq, shift_freq)

        x = freq * Constants.c_in_cm * Constants.h / Constants.k_b
        s += Constants.r * (((x / temp) / (np.exp(x / temp) - 1.0)) -
                                np.log(1.0 - np.exp(-x / temp)))

    return s


def calc_grimme_s_vib(molecule, temp, omega_0, alpha):
    """
    Calculate the entropy according to Grimme's qRRHO method of RR-HO
    interpolation in Chem. Eur. J. 2012, 18, 9955

    :param molecule: (otherm.Molecule)
    :param temp: (float) Temperature in K
    :param omega_0: (float) ω0 parameter
    :param alpha: (float) α parameter
    :return: (float) S_vib in J K-1 mol-1
    """
    s = 0.0

    i_mat = calc_moments_of_inertia(molecule.xyzs)

    # Average I = (I_xx + I_yy + I_zz) / 3.0
    b_avg = np.trace(i_mat) / 3.0

    for freq in molecule.real_vib_freqs():

        omega = freq * Constants.c_in_cm
        mu = Constants.h / (8.0 * np.pi**2 * omega)
        mu_prime = (mu * b_avg) / (mu + b_avg)

        x = freq * Constants.c_in_cm * Constants.h / (Constants.k_b * temp)
        s_v = Constants.r * ((x / (np.exp(x) - 1.0)) - np.log(1.0 - np.exp(-x)))
        s_r = Constants.r * (0.5 + np.log(np.sqrt((8.0 * np.pi**3 * mu_prime * Constants.k_b * temp) /
                                                  (Constants.h**2)
                                                  )))

        w = 1.0 / (1.0 + (omega_0 / freq)**alpha)

        s += w * s_v + (1.0 - w) * s_r

    return s


def calc_entropy(molecule, method='grimme', temp=298.15, ss='1M',
                 shift=100, w0=100, alpha=4):
    """
    Calculate the entropy

    :param molecule: (otherm.Molecule)

    :return: (float) S in J K-1 mol-1
    """

    # Translational entropy component
    s_trans = calc_s_trans_pib(molecule, ss=ss, temp=temp)

    if molecule.n_atoms == 1:
        # A molecule with only one atom has no rotational/vibrational DOF
        return s_trans

    # Rotational entropy component
    s_rot = calc_s_rot_rr(molecule, temp=temp)

    # Vibrational entropy component
    if method.lower() == 'igm':
        s_vib = calc_igm_s_vib(molecule, temp)

    elif method.lower() == 'truhlar':
        s_vib = calc_truhlar_s_vib(molecule, temp, shift_freq=shift)

    elif method.lower() == 'grimme':
        s_vib = calc_grimme_s_vib(molecule, temp, omega_0=w0, alpha=alpha)

    else:
        raise NotImplementedError

    return s_trans + s_rot + s_vib


def calc_zpe(molecule):
    """
    Calculate the zero point energy of a molecule, contributed to by the real
    (positive) frequencies

    :param molecule: (otherm.Molecule)
    :return: (float) E_ZPE
    """

    zpe = 0.0

    for freq in molecule.real_vib_freqs():
        zpe += 0.5 * Constants.h * Constants.n_a * Constants.c_in_cm * freq

    return zpe


def calc_internal_vib_energy(molecule, temp):
    """
    Calculate the internal energy from vibrational motion within the IGM

    :param molecule: (otherm.Molecule)
    :param temp: (float)
    :return: (float) U_vib
    """

    e_vib = 0.0

    # Final 6 vibrational frequencies are translational/rotational
    for freq in molecule.freqs[:-6]:
        if freq > 0:
            x = freq * Constants.c_in_cm * Constants.h / Constants.k_b
            e_vib += Constants.r * x * (1.0 / (np.exp(x/temp) - 1.0))

    return e_vib


def calc_internal_energy(molecule, temp):
    """
    Calculate the internal energy of a molecule

    :param molecule: (otherm.Molecule)
    :param temp: (float) Temperature in K
    :return: (float) U
    """

    zpe = calc_zpe(molecule)
    e_trns = 1.5 * Constants.r * temp

    if is_linear(molecule.coords()):
        # Linear molecules only have two rotational degrees of freedom -> RT
        e_rot = Constants.r * temp

    else:
        # From equipartition with 3 DOF -> 3/2 RT contribution to the energy
        e_rot = 1.5 * Constants.r * temp

    e_vib = calc_internal_vib_energy(molecule, temp=temp)

    return molecule.e + zpe + e_trns + e_rot + e_vib


class Molecule:

    def shift_to_com(self):
        """Shift a molecules xyzs to the center of mass"""

        shifted_xyzs = []

        for xyz_line in self.xyzs:
            pos = np.array(xyz_line[1:]) - self.com
            shifted_xyzs.append(xyz_line[:1] + pos.tolist())

        self.xyzs = shifted_xyzs
        self.com = self.calculate_com()
        return None

    def real_vib_freqs(self, make_real=False):
        """Return the real (positive) vibrational frequencies"""
        # Vibrational frequencies are all but the 6 smallest (rotational +
        # translational) and also remove the largest imaginary frequency if
        # this species is a transtion state
        excluded_n = 7 if self.is_ts else 6

        # Frequencies are sorted high -> low(negative)
        if make_real:
            return [np.abs(freq) for freq in self.freqs[:-excluded_n]]

        else:
            return [freq for freq in self.freqs[:-excluded_n] if freq > 0]

    def calculate_mass(self):
        """Calculate the molecular mass of this molecule in kg"""

        atomic_symbols = [xyz[0] for xyz in self.xyzs]
        masses_amu = [Constants.atomic_masses[elm] for elm in atomic_symbols]

        return Constants.amu_to_kg * sum(masses_amu)

    def calculate_com(self):
        """
        Calculate the center of mass (COM

        :return: (np.ndarray) COM vector
        """
        total_mass = self.mass / Constants.amu_to_kg

        com_vec = np.zeros(3)  # Blank 3D vector for COM vector

        for xyz_line in self.xyzs:
            r_vec = np.array(xyz_line[1:])  # Vector for that atom
            mass = Constants.atomic_masses[xyz_line[0]]
            com_vec += (1.0 / total_mass) * mass * r_vec

        return com_vec

    def coords(self):
        """Return a numpy array shape (n_atoms, 3) of (x,y,z) coordinates"""
        return np.array([np.array(line[1:4]) for line in self.xyzs])

    def calculate_thermochemistry(self, temp=298.15, ss='1M', method='grimme',
                                  shift=100, w0=100, alpha=4, calc_sym=True,
                                  symm_n=None, make_freqs_real=True):
        """
        Calculate thermochemical components and the energies U, H, S, G

        -----------------------------------------------------------------------
        :param temp: (float) Temperature in K

        :param ss: (str) standard state e.g. 1M or 1atm

        :param method: (str) Method to calculate the entropy

        :param shift: (float) Shift frequency used in the Truhlar method of
                      calculating vibrational entropy. All harmonic freqencies
                      below this value will be shifted to this value

        :param w0: (float) ω0 parameter in the Grimme vibrational entropy
                   method

        :param alpha: (float) α parameter the Grimme vibrational entropy
                   method

        :param calc_sym: (bool) Force the calculation of symmetry

        :param symm_n: (int) Override the calculated symmetry number

        :param make_freqs_real: (bool) Convert all imaginary frequencies to
                                their real (positive) analogues
        """
        assert len(self.freqs) == 3 * len(self.xyzs)

        # If the calculation of rotational symmetry number σR is requested or
        # there aren't too many atoms
        if calc_sym or self.n_atoms < 50:
            self.sigma_r = calc_symmetry_number(self)

        # Allow overwriting σR
        if symm_n:
            self.sigma_r = symm_n

        # Convert all imaginary frequencies printed by ORCA as negative to real
        # (so positive) values
        if make_freqs_real:
            self.freqs = [np.abs(freq) for freq in self.freqs]

        self.s = calc_entropy(self, method, temp, ss, shift, w0, alpha)
        self.u = calc_internal_energy(self, temp)
        self.h = self.u + Constants.r * temp
        self.g = self.h - temp * self.s

        return None

    def __init__(self, filename, is_ts=False):
        """
        Molecule initialised from an ORCA output file

        :param filename: (str)
        """
        # Is this molecule a transition state, and so should have one imaginary
        # (negative frequency) that does not contribute to the energy?
        self.is_ts = is_ts

        # Harmonic vibrational frequencies in cm-1
        self.freqs = extract_frequencies(filename)

        # Atom positions [[atom, x, y, z], ...] x/y/z in Å
        self.xyzs = extract_xyzs(filename)
        self.n_atoms = len(self.xyzs)

        # Mass in kg
        self.mass = self.calculate_mass()

        # Matrix of I values in kg m^2
        self.moments_of_inertia = calc_moments_of_inertia(self.xyzs)

        # Centre of mass np.array shape (3,) x/y/z in Å
        self.com = self.calculate_com()
        self.shift_to_com()

        # Rotational symmetry number
        self.sigma_r = 1

        # Total electronic (E), internal (U), enthalpy (H), entropy (S) and
        # Gibbs (free) energy (G) all in molar SI units
        self.e = extract_final_electronic_energy(filename)
        self.u = None
        self.h = None
        self.s = None
        self.g = None


if __name__ == '__main__':

    args = get_args()
    mol = Molecule(args.filename, is_ts=args.transition_state)

    mol.calculate_thermochemistry(temp=args.temp,
                                  ss=args.standard_state,
                                  method=args.method,
                                  shift=args.shift,
                                  w0=args.w0,
                                  alpha=args.alpha,
                                  calc_sym=args.calc_sym,
                                  symm_n=args.symn,
                                  make_freqs_real=args.real_freqs)
    print_output(mol)
