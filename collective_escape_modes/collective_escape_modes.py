# BigBoro - Arbitrary-precision Python/Go software for Böröczky packings with ECMC computations
# https://github.com/jellyfysh/BigBoro
# Copyright (C) 2021 Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, and Werner Krauth
#
# This file is part of BigBoro.
#
# BigBoro is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# BigBoro is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with BigBoro in the LICENSE file.
# If not, see <https://www.gnu.org/licenses/>.
#
# If you use BigBoro in published work, please cite the following reference (see [Hoellmer2021] in References.bib):
# Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, and Werner Krauth,
# Sparse hard-disk packings and Markov chains,
# arXiv e-prints: 2109.13343 (2021), https://arxiv.org/abs/2109.13343.
#
"""Executable script that computes collective escape modes of two-dimensional hard-disk packings."""
import argparse
import sys
from typing import Optional, Sequence
import warnings
import numpy as np


def print_start_message() -> None:
    """"Print the start message which includes the copyright."""
    print("BigBoro (version 1.0) - Arbitrary-precision Python/Go software for Böröczky packings with ECMC "
          "computations - https://github.com/jellyfysh/BigBoro")
    print("Copyright (C) 2021 Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, and Werner Krauth")
    print()


def parse_options(args: Sequence[str]) -> argparse.Namespace:
    """
    Parse command line arguments from the argument strings and store them as attributes of the argparse namespace.

    This function adds the filename, which stores the two-dimensional hard-disk packing in a periodic square box,
    as an obligatory positional argument.

    This function also adds the following (optional) command line arguments:
    1. --length LENGTH  specify the system length of the periodic square box; if None, try to parse the system length
                        from the file that stores the packing (default=None)
    2. --file FILE      set the output filename for the collective escape modes; if empty, the escape modes are stored
                        in 'EscapeModes.txt' (default ='')
    3. --plot PLOT      set the pdf-filename for the plots of all escape modes (requires matplotlib); if empty, no plot
                        is created (default='')
    4. --version        show program's version number and exit
    5. -h, --help       show the help message and exit

    The argument strings can be, for example, sys.argv[1:] in order to parse the command line arguments.

    Parameters
    ----------
    args : Sequence[str]
        The argument strings.

    Returns
    -------
    argparse.Namespace
        The populated argparse namespace.
    """
    parser = argparse.ArgumentParser(
        description="compute collective escape modes of two-dimensional hard-disk packings in a periodic square box")
    parser.add_argument("packing_file", help="specify the path to the file that stores the packing")
    parser.add_argument("--length", type=float, default=None,
                        help="specify the system length of the periodic square box; if None, try to parse the system "
                             "length from the file that stores the packing (default=None)")
    parser.add_argument("--file", default="", help="set the output filename for the collective escape modes; if empty, "
                                                   "the escape modes are stored in 'EscapeModes.txt' (default ='')")
    parser.add_argument("--plot", default="", help="set the pdf-filename for the plots of all escape modes (requires "
                                                   "matplotlib); if empty, no plot is created (default='')")
    parser.add_argument("--version", action="version", version="BigBoro software package version 1.0.")
    args = parser.parse_args(args)
    if args.plot != "" and not args.plot.endswith(".pdf"):
        parser.error(f"argument --plot: invalid value {args.plot}; the filename for the plots of all escape modes has "
                     f"to have the suffix '.pdf'")
    return args


def read_packing(filename: str, system_length: Optional[float]) -> (np.ndarray, float):
    """
    Parse the two-dimensional hard-disk packing in a square periodic box from the given file and return it together with
    the system length of the box.

    The file should have two columns that store the x- and y- positions of each hard disk, respectively. The packing is
    returned as a numpy array of shape (number_of_disks, 2).

    If the system_length argument is None, this function attempts to parse the system length from the file. For this, it
    expects a header line of the format '# System length: {system_length}', where {system_length} is a float. Such a
    line is, for example, included by the construct_packing.py script of the BigBoro package.

    This function returns the (possibly parsed) system length together with the hard-disk packing.

    Parameters
    ----------
    filename : str
        The filename that stores the two-dimensional hard-disk packing.
    system_length : float or None
        The system length of the periodic square box; if None, parse the system length from the packing file.

    Returns
    -------
    (np.ndarray[dtype=np.float64], float)
        The two-dimensional hard-disk packing, the system length.

    Raises
    ------
    AssertionError
        If any non-header line, which does not start with '#', does not contain two floats separated by a whitespace.
        If any position along any axis is not in the interval [0, system_length].
    RuntimeError
        If system_length is None but the system length could not be parsed from the file.
    """
    packing = []
    with open(filename, "r") as file:
        if system_length is None:
            for line in file:
                if line.startswith("# System length: "):
                    system_length = float(line[len("# System length: "):])
                    break
        if system_length is None:
            raise RuntimeError(f"[ERROR] Could not parse system length from the packing stored in {filename}.")
        sphere_index = 0
        for line in file:
            if line.startswith("#"):
                continue
            split_line = line.split()
            assert len(split_line) == 2
            packing.append([np.float64(split_line[0]), np.float64(split_line[1])])
            assert 0.0 <= packing[-1][0] <= system_length
            assert 0.0 <= packing[-1][1] <= system_length
            sphere_index += 1
    return np.array(packing), system_length


def separation_vector(position_one: np.ndarray, position_two: np.ndarray, system_length: float) -> np.ndarray:
    """
    Return the shortest separation vector between the two given positions in a periodic square box of the given system
    length.

    Parameters
    ----------
    position_one : np.ndarray[dtype=np.float64]
        The first position.
    position_two : np.ndarray[dtype=np.float64]
        The second position.
    system_length : float
        The system length of the periodic square box.

    Returns
    -------
    np.ndarray[dtype=np.float64]
        The shortest separation vector (possibly corrected for periodic boundary conditions).

    Raises
    ------
    AssertionError
        If the separation vector could not be corrected for periodic boundary conditions (most probably because the
        given positions did not lie within the square box with each side ranging from 0 to system_length).
    """
    diff = position_one - position_two
    for d in range(diff.shape[0]):
        if diff[d] >= system_length / 2.0:
            diff[d] -= system_length
        elif diff[d] < -system_length / 2.0:
            diff[d] += system_length
        assert -system_length / 2.0 <= diff[d] < system_length / 2.0
    return diff


def check_packing(packing: np.ndarray, system_length: float, distance_precision: float) -> np.ndarray:
    """
    Test that the given two-dimensional hard-disk packing has no overlaps in a periodic square box within the given
    distance precision. Return the contacts within the packing.

    Two hard disks overlap if their distance (that is possibly corrected for periodic boundary conditions) is smaller
    than 2 - distance_precision. Two hard disks are in contact if their distance lies within the interval
    [2 - distance_precision, 2 + distance_precision].

    The two-dimensional hard-disk packing should be an numpy array of shape (number_of_disks, 2) that stores the x- and
    y-position of each hard disk in the periodic square box, respectively. The returned numpy array contains the indices
    of the hard disks in contact in each row, and is thus of shape (number_of_contacts, 2).

    Parameters
    ----------
    packing : np.ndarray[dtype=np.float64]
        The positions of the hard disks in the two-dimensional packing.
    system_length : float
        The system length of the periodic square box of the hard-disk packing.
    distance_precision : float
        The distance precision.

    Returns
    -------
    np.ndarray[dtype=int]
        The contacts in the hard-disk packing.

    Raises
    ------
    RuntimeError
        If there are overlaps between hard disks in the packing.
    """
    number_contacts = 0
    number_overlaps = 0
    contacts = []
    for i in range(packing.shape[0]):
        number_contacts_i = 0
        for j in range(packing.shape[0]):
            if i == j:
                continue
            diff = separation_vector(packing[i, :], packing[j, :], system_length)
            distance = np.sqrt(np.dot(diff, diff))
            if distance < 2.0 - distance_precision:
                print(f"[ERROR] Disks with indices {i} and {j} overlap at distance {distance}.", file=sys.stderr)
                number_overlaps += 1
            if distance <= 2.0 + distance_precision:
                number_contacts_i += 1
                # Avoid double counting.
                if j < i:
                    contacts.append([j, i])
                    number_contacts += 1
    if number_overlaps != 0:
        raise RuntimeError(f"[ERROR] Detected overlaps between disks in the Böröczky packing.")
    return np.array(contacts, dtype=int)


def compute_escape_matrix(packing: np.ndarray, contacts: np.ndarray, system_length: float,
                          distance_precision: float) -> np.ndarray:
    r"""
    Return the escape matrix for the given two-dimensional hard-disk packing with the given contacts.

    The given two-dimensional hard-disk packing should be a numpy array of shape (number_of_disks, 2) that stores the x-
    and y-position of each hard disk in the periodic square box with the given system length, respectively. The given
    contacts numpy array contains the indices of the hard disks in contact in each row, and is thus of shape
    (number_of_contacts, 2). Two hard disks are in contact if their distance lies within the interval
    [2 - distance_precision, 2 + distance_precision].

    The escape matrix is the coefficient matrix of the homogeneous system of linear equations
    (\vec{x}_i - \vec{x}_j) \cdot (\vec{\Delta}_i - \vec{\Delta}_j) = 0 for all disks i and j at the positions \vec{x}_i
    and \vec{x}_j that are in contact. The vector \vec{\Delta}_i is the (infinitesimal) displacement of disk i. The
    system of linear equations expresses that the distances between the hard disks should be unchanged to first order in
    the displacement vectors (see eq. (5) of []). TODO: Update reference!

    The returned escape matrix is of shape (number_of_contacts, 2 * number_of_disks). Each row has four non-zero
    elements corresponding to the contact between hard disks i and j:
    ..., (\vec{x}_i-\vec{x}_j)_x, (\vec{x}_i-\vec{x}_j)_y, ..., (\vec{x}_j-\vec{x}_i)_x, (\vec{x}_j-\vec{x}_i)_y, ...

    Parameters
    ----------
    packing : np.ndarray[dtype=np.float64]
        The positions of the hard disks in the two-dimensional packing.
    contacts : np.ndarray[dtype=int]
        The contacts in the hard-disk packing.
    system_length : float
        The system length of the periodic square box of the hard-disk packing.
    distance_precision : float
        The distance precision.

    Returns
    -------
    np.ndarray[dtype=np.float64]
        The escape matrix.

    Raises
    ------
    AssertionError
        If the contacts numpy array contains a pair of hard disks that is not in contact within the given distance
        precision.
    """
    number_contacts = contacts.shape[0]
    number_disks = packing.shape[0]
    escape_matrix = np.zeros(shape=(number_contacts, 2 * number_disks), dtype=np.float64)
    for contact_i in range(number_contacts):
        i = contacts[contact_i, 0]
        j = contacts[contact_i, 1]
        assert i != j
        diff = separation_vector(packing[i, :], packing[j, :], system_length)
        distance = np.sqrt(np.dot(diff, diff))
        assert 2.0 - distance_precision <= distance <= 2.0 + distance_precision
        escape_matrix[contact_i, 2 * i] = diff[0]
        escape_matrix[contact_i, 2 * i + 1] = diff[1]
        escape_matrix[contact_i, 2 * j] = -diff[0]
        escape_matrix[contact_i, 2 * j + 1] = -diff[1]
    return escape_matrix


def compute_escape_modes(escape_matrix: np.ndarray, number_of_disks: int) -> np.ndarray:
    """
    Return an orthogonal basis of the solution space for the homogeneous system of linear equations parameterized by the
    given escape (coefficient) matrix.

    The coefficient matrix corresponds to the escape matrix that demands that the distance between any two hard disks
    in contact stays constant up to first order in the displacement of every hard disk (see 'compute_escape_matrix'
    function for details). The vectors of the orthogonal basis are thus collective escape modes of the hard-disk
    packing.

    We obtain the orthogonal basis of the solution space for the homogeneous system of linear equations
    escape_matrix times \vec{Delta} = \vec{0} by using a singular value decomposition. The basis vectors of the
    solution space correspond to the singular values that are zero.

    The escape modes are stored in a numpy array of shape (number_of_contacts, 2 * number_of_disks). Each row contains
    the displacement vectors \vec{Delta}_i for every disk i as (\vec{Delta}_1)_x, (\vec{Delta}_1)_y,
    (\vec{Delta}_2)_x, "(\vec{Delta}_2)_y, ....

    This function also tests whether the two trivial global displacement vectors that displace all hard disks along the
    x- and y-axis, respectively, belong to the vector space spanned by the orthogonal basis.

    Parameters
    ----------
    escape_matrix : np.ndarray[dtype=np.float64]
        The escape (coefficient) matrix.
    number_of_disks : int
        The number of disks.

    Returns
    -------
    np.ndarray[dtype=np.float64]
        The orthogonal collective escape modes.

    Raises
    ------
    RuntimeError
        If no escape modes are found.
        If escape_matrix times escape_mode != 0 for any escape mode.
        If any of the two trivial global translation modes is not part of the vector space spanned by the escape modes.
    """
    u, s, v = np.linalg.svd(escape_matrix)
    number_small_singular_values = s[s < 1.0e-13].size
    if number_small_singular_values != 0:
        warnings.warn(f"Computed singular values: {s}.\n"
                      f"Number of singular values that are smaller than 1.0e-13: {number_small_singular_values}.\n"
                      f"This script assumes that these small singular values are zero within floating point precision.",
                      RuntimeWarning)
    # s matrix only contains non-zero singular values.
    number_escape_modes = 2 * number_of_disks - s.shape[0] + number_small_singular_values
    if number_escape_modes == 0:
        raise RuntimeError("[ERROR] Found no escape modes.")
    escape_modes = v[-number_escape_modes:]
    for i in range(number_escape_modes):
        if not np.all(np.abs(escape_matrix @ escape_modes[i]) < 1.0e-13):
            raise RuntimeError("[ERROR] The matrix product escape_matrix times escape mode does not vanish.")
    translation_v1 = np.zeros(2 * number_of_disks, dtype=np.float64)
    translation_v2 = np.zeros(2 * number_of_disks, dtype=np.float64)
    for i in range(number_of_disks):
        translation_v1[2 * i] = 1.0
        translation_v2[2 * i + 1] = 1.0
    v1 = np.zeros(2 * number_of_disks, dtype=np.float64)
    v2 = np.zeros(2 * number_of_disks, dtype=np.float64)
    for i in range(number_escape_modes):
        v1 += np.dot(escape_modes[i], translation_v1) * escape_modes[i]
        v2 += np.dot(escape_modes[i], translation_v2) * escape_modes[i]
    if not np.all(np.abs(translation_v1 - v1) < 1.0e-13):
        raise RuntimeError("[ERROR] The global translation mode along the x-axis is not part of the vector space "
                           "spanned by the escape modes.")
    if not np.all(np.abs(translation_v2 - v2) < 1.0e-13):
        raise RuntimeError("[ERROR] The global translation mode along the y-axis is not part of the vector space "
                           "spanned by the escape modes.")
    return escape_modes


def check_escape_modes(packing: np.ndarray, contacts: np.ndarray, system_length: float,
                       escape_modes: np.ndarray, delta: float) -> None:
    r"""
    Displace the hard disks by their displacement vector in the collective escape modes and check if no overlaps are
    created.

    The two-dimensional hard-disk packing should be an numpy array of shape (number_of_disks, 2) that stores the x- and
    y-position of each hard disk in the periodic square box of the given system length, respectively. The given contacts
    numpy array contains the indices of the hard disks in contact in each row, and is thus of shape
    (number_of_contacts, 2). The escape modes should be a numpy array of shape
    (number_of_contacts, 2 * number_of_disks). Each row contains the displacement vectors \vec{Delta}_i for every hard
    disk i as (\vec{Delta}_1)_x, (\vec{Delta}_1)_y, (\vec{Delta}_2)_x, "(\vec{Delta}_2)_y, ....

    The hard disks are displaced by delta * \vec{Delta}_i. The resulting packing is checked to have no overlaps within
    float precision (i.e., the distance_precision of the 'check_packing' function is set to 1.0e-13). This function also
    checks if the contacts in the packing remain the same up to the precision of delta (distance_precision=delta).

    Parameters
    ----------
    packing : np.ndarray[dtype=np.float64]
        The positions of the hard disks in the two-dimensional packing.
    contacts : np.ndarray[dtype=int]
        The contacts in the hard-disk packing.
    system_length : float
        The system length of the periodic square box of the hard-disk packing.
    escape_modes : np.ndarray[dtype=np.float64]
        The orthogonal collective escape modes.
    delta : float
        The multiplicative prefactor of the displacement vectors in the collective escape modes by which the hard disks
        are displaced.

    Raises
    ------
    RuntimeError
        If there are overlaps between hard disks in the packing after the displacements of a collective escape mode
        (via the 'check_packing' function).
        If the contacts between hard disks were changed (up to the precision of delta) after the displacements of a
        collective escape mode.
    """
    number_of_disks = packing.shape[0]
    number_escape_modes = escape_modes.shape[0]
    for i in range(number_escape_modes):
        copy_packing = np.copy(packing)
        for j in range(number_of_disks):
            copy_packing[j, 0] += delta * escape_modes[i, 2 * j]
            copy_packing[j, 1] += delta * escape_modes[i, 2 * j + 1]
            for d in range(2):
                if copy_packing[j, d] < 0.0:
                    copy_packing[j, d] += system_length
                elif copy_packing[j, d] > system_length:
                    copy_packing[j, d] -= system_length
        # Check if there are any overlaps within float precision.
        # Note that, depending on delta, contacts might disappear with this precision (expected at roughly delta^2).
        check_packing(packing, system_length, 1.0e-13)
        # Check if the number of contacts stays the same when the precision of delta is used.
        new_contacts = check_packing(packing, system_length, delta)
        if not np.all(new_contacts == contacts):
            raise RuntimeError(f"[ERROR] Contacts changed in the hard-disk packing for escape mode {i} (original "
                               f"number of contacts: {contacts.shape[0]}; new number of contacts: "
                               f"{new_contacts.shape[0]}).")


def plot_escape_modes(packing: np.ndarray, system_length: float, escape_modes: np.ndarray,
                      filename: str, factor: float) -> None:
    r"""
    Plot all collective escape modes for the given two-dimensional hard-disk packing in a periodic square box.

    The two-dimensional hard-disk packing should be an numpy array of shape (number_of_disks, 2) that stores the x- and
    y-position of each hard disk in the periodic square box of the given system length, respectively.

    The escape modes should be a numpy array of shape (number_of_contacts, 2 * number_of_disks). Each row contains the
    displacement vectors \vec{Delta}_i for every hard disk i as (\vec{Delta}_1)_x, (\vec{Delta}_1)_y,
    (\vec{Delta}_2)_x, "(\vec{Delta}_2)_y, ....

    This function plots the escape modes as arrows. Here, the displacement vectors are multiplied by the given factor
    so that they are better visible.

    Parameters
    ----------
    packing : np.ndarray[dtype=np.float64]
        The positions of the hard disks in the two-dimensional packing.
    system_length : float
        The system length of the periodic square box of the hard-disk packing.
    escape_modes : np.ndarray[dtype=np.float64]
        The orthogonal collective escape modes.
    filename : str
        The filename of the pdf that is generated for the plot of the escape modes.
    factor : float
        The multiplicative factor for the plotted arrows corresponding to the displacement vectors in the escape modes.
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    number_of_disks = packing.shape[0]
    number_escape_modes = escape_modes.shape[0]
    with PdfPages(filename) as pdf:
        for i in range(number_escape_modes):
            plt.figure()
            ax = plt.gca()
            for j in range(number_of_disks):
                circ = plt.Circle((packing[j, 0], packing[j, 1]), 1, facecolor="none", edgecolor="k")
                ax.add_patch(circ)
                ax.arrow(x=packing[j, 0], y=packing[j, 1],
                         dx=factor * escape_modes[i, 2 * j], dy=factor * escape_modes[i, 2 * j + 1])
            plt.xlim(0, system_length)
            plt.ylim(0, system_length)
            ax.set_aspect("equal")
            plt.title(f"Escape mode {i}")
            pdf.savefig()
            plt.close()


def main() -> None:
    """Compute the collective escape modes of a hard-disk packing based on the command-line arguments."""
    print_start_message()
    args = parse_options(sys.argv[1:])
    packing, system_length = read_packing(args.packing_file, args.length)
    number_of_disks = packing.shape[0]
    print(f"Number of disks in the two-dimensional hard-disk packing: {number_of_disks}")
    print(f"System length of the periodic square box: {system_length}")
    contacts = check_packing(packing, system_length, 1.0e-13)
    if contacts.size == 0:
        raise RuntimeError("[ERROR] Found no contacts in the two-dimensional hard-disk packing.")
    number_of_contacts = contacts.shape[0]
    print(f"Number of contacts in the hard-disk packing: {number_of_contacts}")
    escape_matrix = compute_escape_matrix(packing, contacts, system_length, 1.0e-13)
    escape_modes = compute_escape_modes(escape_matrix, number_of_disks)
    number_escape_modes = escape_modes.shape[0]
    print(f"Number of basis vectors for the collective escape modes: {number_escape_modes}")
    check_escape_modes(packing, contacts, system_length, escape_modes, 1.0e-8)
    print("Actual displacement of each hard disk by its displacement vector multiplied by 1.0e-8 yielded no "
          "overlapping hard disks and no changes of contacts (within a distance precision of 1.0e-8).")
    filename = "EscapeModes.txt" if args.file == "" else args.file
    print(f"Storing the collective escape modes in the file {filename}.")
    np.savetxt(filename, escape_modes,
               header=f"Escape modes for the two-dimensional hard-disk packing in a periodic square box in the "
                      f"file {args.packing_file}.\nEach row contains all displacement vectors \\vec{{Delta}}_i for "
                      f"every disk i as (\\vec{{Delta}}_1)_x, (\\vec{{Delta}}_1)_y, (\\vec{{Delta}}_2)_x, "
                      f"(\\vec{{Delta}}_2)_y, ...")
    if args.plot != "":
        print(f"Plotting the escape modes in the pdf file {args.plot}.")
        plot_escape_modes(packing, system_length, escape_modes, args.plot, 8.0)


if __name__ == '__main__':
    main()
