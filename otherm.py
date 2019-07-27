#!/usr/bin/env python3
import argparse
import numpy as np
from copy import deepcopy


class Constants:
    k_b = 1.38064852E-23                # J K-1
    h = 6.62607004E-34                  # J s
    n_a = 6.022140857E23                # molecules mol-1
    ha_to_j = 4.359744650E-18           # J Ha-1
    atm_to_pa = 101325                  # Pa
    dm_to_m = 0.1                       # m
    amu_to_kg = 1.660539040E-27         # Kg
    r = k_b * n_a                       # J K-1 mol-1
    c = 299792458                       # m s-1
    c_in_cm = c * 100                   # cm s-1
    ang_to_m = 1E-10                    # m
    # TODO complete this list
    atomic_masses = {
        'H': 1.0079, 'He': 4.0026, 'Li': 6.941, 'Be': 9.0122, 'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.994,
        'F': 18.9984, 'Ne': 20.1797, 'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738,
        'S': 32.065, 'Cl': 35.453, 'K': 39.0983, 'Ar': 39.948, 'Ca': 40.078, 'Sc': 44.9559, 'Ti': 47.867,
        'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938, 'Fe': 55.845, 'Ni': 58.6934, 'Co': 58.9332, 'Cu': 63.546,
        'Zn': 65.39, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.8, 'I': 126.9045,
        'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.078, 'Au': 196.9665, 'Pd': 106.42
    }


def get_args():
    """
    Get arguments to the tool with argparse
    :return: The arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", action='store',
                        help='.out file with freq calculation performed')
    parser.add_argument('-t', '--temp', type=float, default=298,
                        help="Temperature (K) at which to calculate G, H and S. Default: %(default)s")
    parser.add_argument('-ss', '--standard_state', type=str, default='1M',
                        help="Standard state. 1atm for gas phase and 1M for solution phase. Default: %(default)s")
    parser.add_argument('-m', '--method', default='grimme', nargs='?', choices=['igm', 'truhlar', 'grimme'],
                        help='Method by which to calculate G, H and S. Default: %(default)s')
    parser.add_argument('-s', '--shift', type=float, default=100,
                        help="Cut-off (in cm-1) to use in Truhlar's method. Frequencies below "
                             "this will be shifted to this value. Default: %(default)s")
    parser.add_argument('-w', '--w0', type=float, default=100,
                        help="Value of w0 to use in Grimme's interpolation method "
                             "Chem. Eur. J. 2012, 18, 9955 eqn. 8. Default: %(default)s")
    parser.add_argument('-a', '--alpha', type=float, default=4,
                        help="Value of alpha to use in Grimme's interpolation method "
                             "Chem. Eur. J. 2012, 18, 9955 eqn. 8. Default: %(default)s")
    parser.add_argument('-cs', '--calc_sym', action='store_true', default=False,
                        help="Force calculation of symmetry number (n^3 algorithm) used for n_atoms < 50."
                             " Default: %(default)s")
    parser.add_argument('-sn', '--symn', type=int, default=None,
                        help="Override the symmetry number calculation. Default: %(default)s")
    parser.add_argument('-r', '--real_freqs', action='store_true', default=False,
                        help='Convert any imaginary frequencies to their real counterparts')

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

    print("{:<50s}{:>33s}".format('Filename', molecule.filename))
    print("{:<50s}{:>33.1f}".format('Temperature (K)', molecule.temp))
    print("{:<50s}{:>33s}".format('Standard state is', '1 M' if molecule.sol else '1 atm'))
    if molecule.real_freqs:
        print("{:<50s}{:>33s}".format('Treat imaginary (negative) frequencies as real', 'True'))
    print("{:<50s}{:>33s}".format('Calculating using the method of',
                                  'IGM' if molecule.igm else 'Truhlar' if molecule.truhlar else 'Grimme'))
    if molecule.grimme:
        print("{:<50s}{:>33s}".format('', 'Chem. Eur. J. 2012, 18, 9955'))
    if molecule.truhlar:
        print("{:<50s}{:>33s}".format('', 'J. Phys. Chem. B, 2011, 115, 14556'))
    print()
    print("{:<50s}{:>33s}".format('Symmetry number (σ)', str(molecule.sigma_r)))
    print("{:<50s}{:>33.2f}".format('Molecular weight (amu)', molecule.mass))
    print()
    if molecule.excess_imag_freqs > 0:
        print('---------------------------------------WARNING--------------------------------------')
        print('                   Found more imaginary frequencies than expected: ', molecule.excess_imag_freqs)
        print('------------------------------------------------------------------------------------\n')
    print("{:<50s}{:>33.2f}".format('Total entropy (J K-1 mol-1)', molecule.s))
    print("{:<50s}{:>33.2f}".format('Total enthalpy (J mol-1)', molecule.h))
    print("{:<50s}{:>33.2f}".format('Total free energy (J mol-1)', molecule.g))
    print()
    print("{:<50s}{:>33}".format('For convenience E, H, G in Hartrees', ''))
    print(molecule.e_elec_ha, molecule.h_ha, molecule.g_ha, sep=',')
    print("----------------------------------------------------------------------------------")

    return 0


def extract_keywords(filename):
    """
    Extract the keywords used for the ORCA job
    :param filename: Name of the ORCA output file
    :return: List of keywords
    """

    keywords = []

    try:
        input_file_block = False

        for line in open(filename, 'r'):
            if 'INPUT FILE' in line:
                input_file_block = True
            if '****END OF INPUT****' in line:
                break
            if '!' in line and input_file_block:
                keywords = line.split()[3:]                 # Line is in the form "|  1> ! Opt "

    except IOError as err:
        exit("I/O error({0}): {1}".format(err.errno, err.strerror))

    keywords = [keyword.lower() for keyword in keywords]    # Lower case all of the keywords

    return keywords


def ensure_freq_calc(keyword_list):
    """
    From a list of input file keywords make sure a frequency calculation has been requested
    :param keyword_list: List of keywords
    """

    if 'freq' not in keyword_list and 'numfreq' not in keyword_list:
        return exit('A frequency calculation has not been requested. Exiting')


def extract_freqs(filename):
    """
    Extract the frequencies from an ORCA output file. The function will first ensure a frequency calculation has
    been performed, then go through the reversed file to find the frequencies (in cm-1). This is done in reverse in case
    multiple hessians have been calculated, the final one is required.
    :param filename: Name of the ORCA output file
    :return: List of frequencies (high to low, in cm-1)
    """

    orca_out_file_lines = []

    try:
        orca_out_file_lines = [line for line in open(filename, 'r')]
    except IOError as err:
        exit("I/O error({0}): {1}".format(err.errno, err.strerror))

    freq_block = False
    freq_list = []

    for line in reversed(orca_out_file_lines):

        if 'NORMAL MODES' in line:
            freq_block = True
        if 'VIBRATIONAL FREQUENCIES' in line:
            break
        if 'cm**-1' in line and freq_block:
            try:
                freq_list.append(float(line.split()[1]))             # Line is in the form "   0:         0.00 cm**-1"
            except TypeError:
                exit("Couldn't extract frequencies. Exiting")

    return freq_list


def extract_final_elec_energy_J_mol(filename):
    """
    Get the final electronic energy from an ORCA output file
    :param filename:
    :return:
    """
    elec_energy = 0.0
    orca_out_lines = open(filename, 'r').readlines()

    for line in orca_out_lines:
        if 'FINAL SINGLE POINT ENERGY' in line:
            elec_energy = Constants.ha_to_j * Constants.n_a * float(line.split()[4])

    return elec_energy


def get_excess_imag_freqs(freq_list, keyword_list):
    """
    Check a list of frequencies to make sure they contain the correct number of imaginary frequencies
    :param freq_list: List of frequencies
    :param keyword_list: List f keywords. Used to determine if a TS was desired
    """

    if len(freq_list) == 0:
        exit('No frequencies were found. Exiting')

    ts_opt = True if 'optts' in keyword_list else False

    n_negative_freqs = len([freq for freq in freq_list if freq < 0.0])
    ideal_n_negative_freqs = 0 if not ts_opt else 1

    return n_negative_freqs - ideal_n_negative_freqs


def extract_mass(filename):
    """
    Extract the total mass (in amu) from a ORCA output file
    :param filename: Name of the ORCA output file
    :return: The mass in AMU
    """

    try:
        for line in open(filename, 'r'):
            if 'Total Mass' in line:
                return float(line.split()[-2])
    except IOError as err:
        exit("I/O error({0}): {1}".format(err.errno, err.strerror))


def calc_com(xyz_list):
    """
    From a list of xyzs compute the center of mass
    :param xyz_list: List of xyzs in the format [[C, 0.0000, 0.0000, 0.0000], ....]
    :return: The COM vector as a np array
    """

    atom_names = [line[0] for line in xyz_list]  # List of atom labels [C, H, N...]
    total_mass = sum([Constants.atomic_masses[atom] for atom in atom_names])  # Sum all the atomic masses au
    com_vec = np.zeros(3)  # Blank 3D vector for COM vector

    for xyz_line in xyz_list:
        r_vec = np.array(xyz_line[1:])  # Vector for that atom
        mass = Constants.atomic_masses[xyz_line[0]]
        com_vec += (1.0 / total_mass) * mass * r_vec

    return com_vec


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

    def get_xyzs():

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

    def shift_xyz_to_com(xyz_list):

        com_vec = calc_com(xyz_list)
        shifted_xyzs = []

        for xyz_line in xyz_list:

            atom_label = xyz_line[0]
            r_vec = np.array(xyz_line[1:])

            shifted_xyzs.append([atom_label] + np.ndarray.tolist(r_vec - com_vec))

        return shifted_xyzs

    xyzs = shift_xyz_to_com(get_xyzs())

    return xyzs


def xyz2coord(xyzs):
    """
    For a set of xyzs in the form e.g [[C, 0.0, 0.0, 0.0], ...] convert to a np array of coordinates, containing just
    just the x, y, z coordinates
    :param xyzs: List of xyzs
    :return: numpy array of coords
    """
    if isinstance(xyzs[0], list):
        return np.array([np.array(line[1:4]) for line in xyzs])
    else:
        return np.array(xyzs[1:4])


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


def apply_n_fold_rotation(coords, origin, axis, n):
    """
    This function will SHIFT the origin of the coordinates to (0.0, 0.0, 0.0)

    :param coords:
    :param origin:
    :param axis:
    :param n:
    :return:
    """

    tmp_coords = deepcopy(coords)
    tmp_coords = np.array([coord - origin for coord in tmp_coords])

    rot_mat = rotation_matrix(axis, theta=(2.0 * np.pi / n))
    return np.array([np.matmul(rot_mat, coord) for coord in tmp_coords])


def find_highest_cn(coords, xyzs, origin, max_n, dist_tol):
    """
    Find the highest symmetry rotation axis

    :param coords:
    :param xyzs:
    :param max_n:
    :param dist_tol:
    :return:
    """

    axes = get_possible_axes(coords)
    min_n = 2                                                                           # Minimum n-fold rotation is 2
    n_origin_axis = [[[], []] for _ in range(max_n+1)]                                  # not elegant..

    for axis in axes:
        for n in range(min_n, max_n+1):
            tmp_coords = np.array([coord - origin for coord in deepcopy(coords)])
            tmp_coords_rot = apply_n_fold_rotation(coords, origin, axis, n)

            if is_same_structure(tmp_coords, tmp_coords_rot, xyzs, dist_tol):
                n_origin_axis[n][0].append(origin)
                n_origin_axis[n][1].append(axis)

    max_n, origins_max_n, axes_max_n = 0, [], []
    for n in range(len(n_origin_axis)):
        if len(n_origin_axis[n][0]) > 0:
            origins_max_n, axes_max_n, max_n = n_origin_axis[n][0], n_origin_axis[n][1], n

    return max_n, origins_max_n, axes_max_n


def calc_symmetry_number(coords, xyzs, origin, max_n_fold_rot_searched=6, dist_tol=0.3):
    """
    Caclulate the symmetry number of a molecule. Based on Theor Chem Account (2007) 118:813–826

    :param coords:
    :param xyzs:
    :param origin:
    :param max_n_fold_rot_searched:
    :param dist_tol:
    :return:
    """

    symmetry_number = 1
    max_n, origins, axes = find_highest_cn(coords, xyzs, origin, max_n_fold_rot_searched, dist_tol)

    if max_n == 1:
        return symmetry_number

    possible_axes = get_possible_axes(coords)
    for origin in origins:
        curr_symmetry_number = 1                            # Already has E symmetry
        for axis in possible_axes:
            for n in range(2, max_n_fold_rot_searched+1):

                tmp_coords = np.array([coord - origin for coord in deepcopy(coords)])
                tmp_coords_rot = apply_n_fold_rotation(coords, origin, axis, n)

                if is_same_structure(tmp_coords, tmp_coords_rot, xyzs, dist_tol):
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


def calc_q_trans_igm(molecule):
    """
    Calculate the translational partition function using the PIB model, coupled with an effective volume
    :param molecule: A Molecule object
    """

    effective_volume = 0            # Effective volume in m3

    if molecule.gas:
        effective_volume = Constants.k_b * molecule.temp / Constants.atm_to_pa
    if molecule.sol:
        effective_volume = 1.0 / (Constants.n_a * (1.0 / Constants.dm_to_m)**3)

    molecule.q_trans = ((2.0 * np.pi * molecule.mass_kg * Constants.k_b * molecule.temp / Constants.h**2)**1.5 *
                        effective_volume)

    return molecule.q_trans


def calc_q_rot_igm(molecule):
    """
    Calculate the rotational partition function using the IGM method. Uses the rotational symmetry number, default = 1
    :param molecule: A Molecule object
    """

    omega_x = Constants.h**2 / (8.0 * np.pi**2 * Constants.k_b * molecule.moments_of_inertia[0, 0])
    omega_y = Constants.h**2 / (8.0 * np.pi**2 * Constants.k_b * molecule.moments_of_inertia[1, 1])
    omega_z = Constants.h**2 / (8.0 * np.pi**2 * Constants.k_b * molecule.moments_of_inertia[2, 2])

    if molecule.n_atoms == 1:
        molecule.q_rot = 1
        return molecule.q_rot
    else:
        molecule.q_rot = (np.sqrt(np.pi)/molecule.sigma_r) * (molecule.temp**1.5 / np.sqrt(omega_x * omega_y * omega_z))

    return molecule.q_rot


def calc_q_vib_igm(molecule):
    """
    Calculate the vibrational partition function using the IGM method. Uses the rotational symmetry number, default = 1
    :param molecule: A Molecule object
    """

    temp = molecule.temp

    if molecule.n_atoms == 1:
        molecule.q_vib = 1
        return molecule.q_vib

    for freq in molecule.freqs:
        if freq > 0:                # Imaginary frequencies are discarded
            x = freq * Constants.c_in_cm * Constants.h / Constants.k_b
            molecule.q_vib *= np.exp(-x / (2.0 * temp)) / (1.0 - np.exp(-x / temp))

    return molecule.q_vib


def calc_s_trans_pib(molecule):

    q_trans = calc_q_trans_igm(molecule)
    return Constants.r * (np.log(q_trans) + 1.0 + 1.5)


def calc_s_rot_rr(molecule):

    if molecule.n_atoms == 1:
        return 0

    q_rot = calc_q_rot_igm(molecule)
    return Constants.r * (np.log(q_rot) + 1.5)


def calc_igm_entropy(molecule):
    """
    Calculate the entropy of a molecule according to the Ideal Gas Model (IGM) method
    :param molecule: A Molecule object
    :return: The sum of translational, rotational and vibrational entropies
    """

    s_trans = calc_s_trans_pib(molecule)

    if molecule.n_atoms == 1:
        return s_trans

    s_rot = calc_s_rot_rr(molecule)

    vib_freqs = molecule.freqs[:-6]                       # Vibrational frequencies are all but the 6 lowest eigenvalues
    s_vib = 0.0                                           # which correspond to translation and rotation (3+3)

    for freq in vib_freqs:
        if freq > 0:                # Imaginary frequencies don't contribute to the partition function
            x = freq * Constants.c_in_cm * Constants.h / (Constants.k_b * molecule.temp)
            s_vib += Constants.r * ((x / (np.exp(x) - 1.0)) - np.log(1.0 - np.exp(-x)))

    return s_trans + s_rot + s_vib


def calc_truhlar_entropy(molecule):
    """
    Calculate the entropy of a molecule according to the Truhlar's method of shifting low frequency modes
    :param molecule: A Molecule object
    :return: The sum of translational, rotational and vibrational entropies
    """

    s_trans = calc_s_trans_pib(molecule)
    s_rot = calc_s_rot_rr(molecule)

    vib_freqs = molecule.freqs[:-6]                       # Vibrational frequencies are all but the 6 lowest eigenvalues
    s_vib = 0                                             # which correspond to translation and rotation (3+3)
    temp = molecule.temp
    for freq in vib_freqs:

        if freq < molecule.shift_freq:
            freq = molecule.shift_freq

        x = freq * Constants.c_in_cm * Constants.h / Constants.k_b
        s_vib += Constants.r * (((x / temp) / (np.exp(x / temp) - 1.0)) - np.log(1.0 - np.exp(-x / temp)))

    return s_trans + s_rot + s_vib


def calc_grimme_entropy(molecule):
    """
    Calculate the entropy according to Grimme's qRRHO method of RR-HO interpolation in Chem. Eur. J. 2012, 18, 9955
    :param molecule: A Molecule object
    :return: The entropy in J K-1 mol-1
    """

    def calc_s_vib():

        vib_freqs = molecule.freqs[:-6]
        s = 0.0
        b_avg = np.trace(molecule.moments_of_inertia) / 3.0                 # Average I = (I_xx + I_yy + I_zz) / 3.0

        for freq in vib_freqs:

            if freq > 0:

                omega = freq * Constants.c_in_cm
                mu = Constants.h / (8.0 * np.pi**2 * omega)
                mu_prime = (mu * b_avg) / (mu + b_avg)

                x = freq * Constants.c_in_cm * Constants.h / (Constants.k_b * molecule.temp)
                s_v = Constants.r * ((x / (np.exp(x) - 1.0)) - np.log(1.0 - np.exp(-x)))
                s_r = Constants.r * (0.5 + np.log(np.sqrt((8.0 * np.pi**3 * mu_prime * Constants.k_b * molecule.temp) /
                                                          (Constants.h**2)
                                                          )))

                w = 1.0 / (1.0 + (molecule.omega_0 / freq)**molecule.alpha)

                s += w * s_v + (1.0 - w) * s_r

        return s

    s_trans = calc_s_trans_pib(molecule)
    s_rot = calc_s_rot_rr(molecule)
    s_vib = calc_s_vib()

    return s_trans + s_rot + s_vib


def calc_entropy(molecule):
    """
    Calculate the entropy
    :param molecule: A Molecule object
    :return: The entropy in J K-1 mol-1
    """

    s = 0.0

    if molecule.igm:
        s = calc_igm_entropy(molecule)

    if molecule.truhlar:
        s = calc_truhlar_entropy(molecule)

    if molecule.grimme:
        s = calc_grimme_entropy(molecule)

    return s


def calc_zpe(molecule):
    """
    Calculate the zero point energy of a molecule, contributed to by the real (positive) frequencies

    :param molecule:
    :return:
    """

    zpe = 0.0

    for freq in molecule.freqs[:-6]:            # Vibrational frequencies. First 6 are translational/rotational
        if freq > 0:
            zpe += 0.5 * Constants.h * Constants.n_a * Constants.c_in_cm * freq

    return zpe


def calc_internal_vib_energy(molecule):
    """
    Calculate the internal energy from vibrational motion within the IGM

    :param molecule:
    :return:
    """

    e_vib = 0.0
    for freq in molecule.freqs[:-6]:  # Vibrational frequencies. First 6 are translational/rotational
        if freq > 0:
            x = freq * Constants.c_in_cm * Constants.h / Constants.k_b
            e_vib += Constants.r * x * (1.0 / (np.exp(x/molecule.temp) - 1.0))

    return e_vib


def calc_internal_energy(molecule):

    elec = extract_final_elec_energy_J_mol(molecule.filename)
    zpe = calc_zpe(molecule)
    e_trns = 1.5 * Constants.r * molecule.temp                                                             # Translation
    e_rot = Constants.r * molecule.temp if is_linear(molecule.coords) else 1.5 * Constants.r * molecule.temp  # Rotation
    e_vib = calc_internal_vib_energy(molecule)

    return elec + zpe + e_trns + e_rot + e_vib


class Molecule(object):

    def __init__(self, filename, temp, ss, method, shift, w0, alpha, calc_sym, symmn, real_freqs):

        self.filename = filename
        self.name = self.filename[:-4]                                      # Strip the .out extension from the filename
        self.temp = temp                                                    # K
        self.real_freqs = real_freqs

        self.sol = True if ss == '1M' else False
        self.gas = True if ss == '1atm' else False

        self.igm = True if method == 'igm' else False
        self.truhlar = True if method == 'truhlar' else False
        self.grimme = True if method == 'grimme' else False

        self.shift_freq = shift                                             # cm-1
        self.omega_0 = w0                                                   # cm-1
        self.alpha = alpha

        self.keywords = extract_keywords(self.filename)
        ensure_freq_calc(self.keywords)
        self.freqs = extract_freqs(self.filename)                           # cm-1
        self.excess_imag_freqs = get_excess_imag_freqs(self.freqs, self.keywords)
        if real_freqs:
            self.freqs = [np.abs(freq) for freq in self.freqs]

        self.mass = extract_mass(self.filename)                             # AMU
        self.mass_kg = self.mass * Constants.amu_to_kg                      # kg

        self.xyzs = extract_xyzs(self.filename)                             # [[atom, x, y, z], ....] x/y/z in Å
        self.coords = xyz2coord(self.xyzs)
        self.n_atoms = len(self.xyzs)

        self.moments_of_inertia = calc_moments_of_inertia(self.xyzs)        # Matrix of I values in kg m^2
        self.com = calc_com(self.xyzs)                                      # Centre of mass np.array

        self.sigma_r = 1

        if calc_sym or self.n_atoms < 50:
            self.sigma_r = calc_symmetry_number(self.coords, self.xyzs, self.com)
        if symmn:
            self.sigma_r = symmn

        self.q_trans = None
        self.q_rot = None
        self.q_vib = None

        self.s = calc_entropy(self)                 # J K-1 mol-1
        self.u = calc_internal_energy(self)         # J mol-1
        self.h = self.u + Constants.r * temp        # J mol-1
        self.g = self.h - temp * self.s             # J mol-1

        self.e_elec_ha = extract_final_elec_energy_J_mol(filename) / (Constants.ha_to_j * Constants.n_a)
        self.h_ha = self.h / (Constants.ha_to_j * Constants.n_a)
        self.g_ha = self.g / (Constants.ha_to_j * Constants.n_a)


if __name__ == '__main__':

    args = get_args()
    mol = Molecule(args.filename, args.temp, args.standard_state, args.method,
                   args.shift, args.w0, args.alpha, args.calc_sym, args.symn, args.real_freqs)
    print_output(mol)
