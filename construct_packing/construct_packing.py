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
"""Executable script that generates arbitrary precision Böröczky packings."""
import argparse
from decimal import Decimal, getcontext
from enum import auto, Enum
import sys
import warnings
from typing import Callable, Sequence
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

    This function adds the following (optional) command line arguments:
    1. --core {kahle,boro}   specify the core (default='kahle')
    2. --chain {geometric,circle}
                             specify the shape of the convex polygonal chain (default='geometric')
    3. --layers LAYERS       set the number of layers > 0 in the (half-)branches (default=5)
    4. --phi PHI             set the attenuation parameter phi in (0, 1) for the construction of the geometric convex
                             polygonal chain (default='0.7')
    5. --file FILE           set the output filename for the hard-disk packing; if empty, the packing is stored in
                             {core}_Layers{l}.txt, where {core} is the used core, and {l} is the number of layers
                             (default='')
    6. --precision PRECISION
                             set the precision in places that is used for the construction (default=50)
    7. --bisection BISECTION
                             set the precision {b} of the bisection search for g_2^<, where the result will be
                             correct up to 1.0e-{b}; if {b} <= 0, {b} = {p} - 3 is used, with {p} as the precision
                             in places (default=0)
    8. --plot                plot and show the constructed hard-disk packing (requires matplotlib)
    9. --version             show program's version number and exit
    10. -h, --help           show the help message and exit

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
        description="generate two-dimensional locally stable Böröczky packings in a periodic square box")
    parser.add_argument("--core", choices=["kahle", "boro"], default="kahle", help="specify the core (default='kahle')")
    parser.add_argument("--chain", choices=["geometric", "circle"], default="geometric",
                        help="specify the shape of the convex polygonal chain (default='geometric')")
    parser.add_argument("--layers", type=int, default=5,
                        help="set the number of layers > 0 in the (half-)branches (default=5)")
    parser.add_argument("--phi", default="0.7",
                        help="set the attenuation parameter phi in (0, 1) for the construction of the geometric convex "
                             "polygonal chain (default='0.7')")
    parser.add_argument("--file", default="",
                        help="set the output filename for the hard-disk packing; "
                             "if empty, the packing is stored in {core}_Layers{l}.txt, where {core} is the used "
                             "core, and {l} is the number of layers (default='')")
    parser.add_argument("--precision", type=int, default=50,
                        help="set the precision in places that is used for the construction (default=50)")
    parser.add_argument("--bisection", type=int, default=0,
                        help="set the precision {b} of the bisection search for g_2^<, where the result will be "
                             "correct up to 1.0e-{b}; if {b} <= 0, {b} = {p} - 3 is used, with {p} as the precision "
                             "in places (default=0)")
    parser.add_argument("--plot", action="store_true",
                        help="plot and show the constructed hard-disk packing (requires matplotlib)")
    parser.add_argument("--version", action="version", version="BigBoro software package version 1.0.")
    args = parser.parse_args(args)
    if not args.layers > 0:
        parser.error(f"argument --layers: invalid value {args.layers}; the number of layers has to be bigger than zero")
    if not args.precision > 0:
        parser.error(f"argument --precision: invalid value {args.precision}; the precision in places has to be bigger "
                     f"than zero")
    if args.chain == "circle" and args.phi != "0.7":
        parser.error(f"argument --phi: phi should only be set if a geometric convex polygonal chain is specified")
    try:
        if not 0.0 < float(args.phi) < 1.0:
            parser.error(f"argument --phi: invalid value {args.phi}; "
                         f"the value has to be a float between 0.0 and 1.0")
    except ValueError:
        parser.error(f"argument --phi: invalid float value: '{args.phi}'")
    return args


def pi() -> Decimal:
    """
    Compute pi to the current precision.

    See https://docs.python.org/3/library/decimal.html#recipes.

    Returns
    -------
    decimal.Decimal
        Pi to the current precision.
    """
    getcontext().prec += 2  # extra digits for intermediate steps
    lasts, t, s, n, na, d, da = 0, Decimal(3), 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n + na, na + 8
        d, da = d + da, da + 32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return Decimal(+s)  # unary plus applies the new precision


def cos(x: Decimal) -> Decimal:
    """
    Return the cosine of the angle x as measured in radians.

    The Taylor series approximation works best for a small value of x. For larger values, first compute
    x = x % (2 * pi). See https://docs.python.org/3/library/decimal.html#recipes.

    Parameters
    ----------
    x : Decimal
        The angle in radians.

    Returns
    -------
    Decimal
        The cosine of x.
    """
    x = x % (Decimal(2) * pi())
    getcontext().prec += 2  # extra digits for intermediate steps
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return Decimal(+s)  # unary plus applies the new precision


def sin(x: Decimal) -> Decimal:
    """
    Return the sine of the angle x as measured in radians.

    The Taylor series approximation works best for a small value of x. For larger values, first compute
    x = x % (2 * pi). See https://docs.python.org/3/library/decimal.html#recipes.

    Parameters
    ----------
    x : Decimal
        The angle in radians.

    Returns
    -------
    Decimal
        The sine of x.
    """
    x = x % (Decimal(2) * pi())
    getcontext().prec += 2  # extra digits for intermediate steps
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i - 1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return Decimal(+s)  # unary plus applies the new precision


def asin(x: Decimal) -> Decimal:
    """
    Return the arcsine as measured in radians of of x.

    Parameters
    ----------
    x : Decimal
        The argument where -1 <= x <= 1.

    Returns
    -------
    Decimal
        The arcsine of x.

    Raises
    ------
    AssertionError
        If not -1 <= x <= 1.
    """
    assert -1 <= x <= 1
    if abs(x) < Decimal(1) / Decimal(2).sqrt():
        getcontext().prec += 2
        xc = x * x
        k = x
        s1 = k
        n = 0
        while True:
            n += 1
            k *= xc * (2 * n - 1) * (2 * n - 1) / (2 * n * (2 * n + 1))
            s2 = s1 + k
            if s2 == s1:
                break
            s1 = s2
        getcontext().prec -= 2
        return s2
    else:
        if x >= 0:
            return pi() / Decimal(2) - asin((Decimal(1) - x * x).sqrt())
        else:
            return -pi() / Decimal(2) + asin((Decimal(1) - x * x).sqrt())


class Core(Enum):
    """
    Class that stores the supported cores of the Böröczky packings as an enumeration.

    Possible cores are BORO and KAHLE.
    """
    BORO = auto()
    KAHLE = auto()


def place_core(core: Core) -> np.ndarray:
    """
    Place the specified core.

    There are 20 disks in the Böröczky core so that the shape of the returned numpy array is (20, 2). In contrast,
    there are 8 disks in the Kahle core so that the shape of the returned numpy array is (8, 2).

    Parameters
    ----------
    core : Core
        The core that should be placed.

    Returns
    -------
    np.ndarray[dtype=Decimal]
        The positions of the disks
    """
    if core is Core.BORO:
        pts_core = []
        sqrt_two = Decimal(2).sqrt()
        two_plus_sqrt_two = Decimal(2) + sqrt_two
        diagonal_value_minus_two = (Decimal(6).sqrt() + Decimal(2).sqrt()) / Decimal(2)
        diagonal_value = diagonal_value_minus_two + Decimal(2)

        pts_core.append([diagonal_value, diagonal_value_minus_two])
        pts_core.append([-diagonal_value, diagonal_value_minus_two])
        pts_core.append([diagonal_value, -diagonal_value_minus_two])
        pts_core.append([-diagonal_value, -diagonal_value_minus_two])
        pts_core.append([diagonal_value_minus_two, diagonal_value])
        pts_core.append([-diagonal_value_minus_two, diagonal_value])
        pts_core.append([diagonal_value_minus_two, -diagonal_value])
        pts_core.append([-diagonal_value_minus_two, -diagonal_value])

        pts_core.append([diagonal_value, diagonal_value])
        pts_core.append([diagonal_value, -diagonal_value])
        pts_core.append([-diagonal_value, -diagonal_value])
        pts_core.append([-diagonal_value, diagonal_value])

        pts_core.append([0, sqrt_two])
        pts_core.append([0, -sqrt_two])
        pts_core.append([-sqrt_two, 0])
        pts_core.append([sqrt_two, 0])

        pts_core.append([0, two_plus_sqrt_two])
        pts_core.append([0, -two_plus_sqrt_two])
        pts_core.append([-two_plus_sqrt_two, 0])
        pts_core.append([two_plus_sqrt_two, 0])
    else:
        assert core is Core.KAHLE
        pts_core = []
        one = Decimal(1)
        one_plus_sqrt_three = Decimal(3).sqrt() + Decimal(1)

        pts_core.append([one, one])
        pts_core.append([one, -one])
        pts_core.append([-one, -one])
        pts_core.append([-one, one])

        pts_core.append([0, one_plus_sqrt_three])
        pts_core.append([0, -one_plus_sqrt_three])
        pts_core.append([-one_plus_sqrt_three, 0])
        pts_core.append([one_plus_sqrt_three, 0])
    return np.array(pts_core)


class Chain(Enum):
    """
    Class that stores the supported convex polygonal chains of the A disks in the Böröczky packings as an enumeration.

    Possible chains are GEOMETRIC and CIRCLE.
    """
    GEOMETRIC = auto()
    CIRCLE = auto()


def place_a_geometric(k_layers: int, delta: Decimal, core: Core, phi: Decimal) -> np.ndarray:
    r"""
    Place k_layers disks on a convex polygonal chain between g_3 = 2 + \sqrt{3} and g_2^< = 2 \sqrt{3} - delta with
    distances 2 between each disk using a geometric progression.

    This function constructs the A disks in the (half-)branch up to the boundary of the system box for which its
    symmetry axis along the branch coincides with the positive x-axis.

    The y-component of the first A_1 disk is g_3 = 2 + \sqrt{3}, while the x-component depends on the used core. Every
    subsequent y-component of disk A_{i+1} is computed via the geometric progression
    dist(A_{i+1}, g_2^<) = phi * dist(A_i, g_2^<), where dist(A_i, g_2^<) is the (shortest) distance between disk A_i
    and the line at height g_2^<.

    The returned numpy array contains all two-dimensional positions of the k_layers disks and is of shape (k_layers, 2).

    Parameters
    ----------
    k_layers : int
        The number of A spheres on the convex polygonal chain.
    delta : Decimal
        The distance that the line g_2^< lies below g_2 = 2 \sqrt{3}.
    core : Core
        The core of the Böröczky packing.
    phi : Decimal
        The attenuation parameter of the geometric progression.

    Returns
    -------
    np.ndarray[dtype=Decimal]
        The positions of the A disks.
    """
    g_2 = Decimal(2) * Decimal(3).sqrt()
    g_3 = Decimal(2) + Decimal(3).sqrt()
    pts_a = np.empty([k_layers, 2], dtype=np.dtype(Decimal))
    pts_a[0, 1] = g_3
    if core is Core.BORO:
        pos_core_disk = Decimal(2) + (Decimal(6).sqrt() + Decimal(2).sqrt()) / Decimal(2)
        pts_a[0, 0] = pos_core_disk + (Decimal(4) - (g_3 - pos_core_disk) ** Decimal(2)).sqrt()
    else:
        assert core is Core.KAHLE
        pts_a[0, 0] = g_3
    for i in range(1, k_layers):
        diff_y = pts_a[i - 1, 1] - g_2 + delta
        angle = asin(Decimal(0.5) * diff_y * (phi - Decimal(1)))
        pts_a[i, 0] = pts_a[i - 1, 0] + Decimal(2) * cos(angle)
        pts_a[i, 1] = pts_a[i - 1, 1] + Decimal(2) * sin(angle)
    return pts_a


def place_a_circle(k_layers: int, delta: Decimal, core: Core) -> np.ndarray:
    r"""
    Place k_layers disks on a convex polygonal chain between g_3 = 2 + \sqrt{3} and g_2^< = 2 \sqrt{3} - delta with
    distances 2 between each disk using a circle.

    This function constructs the A disks in the upper half plane in the (half-)branch up to the boundary of the system
    box for which its symmetry axis along the branch coincides with the positive x-axis.

    The y-component of the first A_1 disk is g_3 = 2 + \sqrt{3}, while the x-component depends on the used core.
    The k_layers disks lie on a circle with a center at distance 1 right to the last disk A_{k_layers} (i.e., on
    the system boundary) with a given radius that, in our case, directly depends on the distance delta that the line
    g_2^< lies below g_2 = 2 \sqrt{3}.

    The returned numpy array contains all two-dimensional positions of the k_layers disks and is of shape (k_layers, 2).

    Parameters
    ----------
    k_layers : int
        The number of A spheres on the convex polygonal chain.
    delta : Decimal
        The distance that the line g_2^< lies below g_2 = 2 \sqrt{3}.
    core : Core
        The core of the Böröczky packing.

    Returns
    -------
    np.ndarray[dtype=Decimal]
        The positions of the A disks.
    """
    g_2 = Decimal(2) * Decimal(3).sqrt()
    g_3 = Decimal(2) + Decimal(3).sqrt()
    half = Decimal('0.5')
    pts_a = np.empty([k_layers, 2], dtype=np.dtype(Decimal))
    pts_a[0, 1] = g_3
    if core is Core.BORO:
        pos_core_disk = Decimal(2) + (Decimal(6).sqrt() + Decimal(2).sqrt()) / Decimal(2)
        pts_a[0, 0] = pos_core_disk + (Decimal(4) - (g_3 - pos_core_disk) ** Decimal(2)).sqrt()
    else:
        assert core is Core.KAHLE
        pts_a[0, 0] = g_3
    diff = g_3 - g_2 + delta
    # The scaling of the (large) radius is approximately length_branch^2 / (2.0 * diff), where length_branch and diff
    # are the distances of the x- and y-components of A_1 and B_{k_layers}, respectively. Here, we approximate
    # length_branch by 2 * k_layers, and assume that g_2^< is tangent to the circle at the lowest point on the circle.
    radius = Decimal(2) * Decimal(k_layers ** 2) / diff
    dtheta = Decimal(2) * asin(Decimal(1) / radius)
    theta = Decimal(k_layers - 1) * dtheta + half * dtheta
    length_branch = radius * sin(theta)
    pts_a[-1, 0] = length_branch + pts_a[0, 0] - 1
    pts_a[-1, 1] = g_3 + radius * cos(theta) - radius * cos(half * dtheta)
    for i in range(2, k_layers):
        pts_a[-i, 0] = pts_a[-i + 1, 0] - Decimal(2) * cos(Decimal(i - 1) * dtheta)
        pts_a[-i, 1] = pts_a[-i + 1, 1] + Decimal(2) * sin(Decimal(i - 1) * dtheta)
    return pts_a


def place_b(a_x: Decimal, a_y: Decimal, c_x: Decimal) -> (Decimal, Decimal):
    """
    Place disk B_i so that is has distance 2 to disk A_i and distance two to disk C_{i-1}.

    This function constructs the B_i disk in the (half-)branch up to the boundary of the system box for which its
    symmetry axis along the branch coincides with the positive x-axis. In this branch, the y-components of the C disks
    are all zero in this branch.

    Parameters
    ----------
    a_x : Decimal
        The x-component of disk A_i.
    a_y : Decimal
        The y-component of disk A_i.
    c_x : Decimal
        The x-component of disk C_{i-1}.

    Returns
    -------
    (Decimal, Decimal)
        The x- and y-component of disk B_i.
    """
    a = (c_x - a_x) / a_y
    b = (a_x ** Decimal(2) + a_y ** Decimal(2) - c_x ** Decimal(2)) / (Decimal(2) * a_y)
    delta_div4 = (Decimal(4) - Decimal(2) * a * b * c_x + Decimal(4) * a ** Decimal(2)
                  - b ** Decimal(2) - c_x ** Decimal(2) * a ** Decimal(2))
    b_x = (c_x - a * b + delta_div4.sqrt()) / (Decimal(1) + a ** Decimal(2))
    b_y = a * b_x + b
    return b_x, b_y


def place_c(b_x: Decimal, b_y: Decimal) -> Decimal:
    """
    Place disk C_i so that is has distance 2 to disk B_i.

    This function constructs the C_i disk in the (half-)branch up to the boundary of the system box for which its
    symmetry axis along the branch coincides with the positive x-axis. In this branch, the y-components of the C disks
    are all zero.

    Parameters
    ----------
    b_x : Decimal
        The x-component of disk B_i.
    b_y : Decimal
        The y-component of disk B_i.

    Returns
    -------
    Decimal
        The x-component of disk C_i.
    """
    delta_div4 = Decimal(4) - b_y ** Decimal(2)
    c_x = b_x + delta_div4.sqrt()
    return c_x


def place_bc(k_layers: int, pts_a: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Place the B and C disks in the (half-)branch with k_layers layers based on the A disks on a convex polygonal
    chain.

    This function constructs the B and C disks in the upper half plane in the (half-)branch up to the boundary of the
    system box for which its symmetry axis along the branch coincides with the positive x-axis.

    The disk B_1 lies below disk A_1. The function then iteratively places C_i on the x-axis to the right of B_i at
    distance 2, and B_{i+1} to the right of A_i so that it has distance 2 to the disks A_{i+1} and C_i.

    The function places k_layers B disks but only k_layers - 1 C disks in order to allow for periodic boundary
    conditions with an vertical additional symmetry axis at x = L / 2 going through B_{k_layers}. (This only succeeds,
    however, if the A disks were placed with the correct line g_2^< so that the distance in the x-component between
    A_{k_layers} and B_{k_layers} is 1. This is not checked in this function.)

    The returned numpy arrays contain all two-dimensional positions of the B and C disks, respectively, and are of
    shape (k_layers, 2) and (k_layers - 1, 2).

    Parameters
    ----------
    k_layers : int
        The number of B spheres (and the number of C spheres + 1).
    pts_a : np.ndarray[dtype=Decimal]
        The positions of the A disks in a numpy array of shape (k_layers, 2).

    Returns
    -------
    (np.ndarray[dtype=Decimal], np.ndarray[dtype=Decimal])
        The positions of the B disks, the positions of the C disks.
    """
    pts_b = np.empty([k_layers, 2], dtype=np.dtype(Decimal))
    pts_c = np.empty([k_layers - 1, 2], dtype=np.dtype(Decimal))
    pts_b[0, 0] = pts_a[0, 0]
    pts_b[0, 1] = pts_a[0, 1] - Decimal(2)
    for i in range(1, k_layers):
        pts_c[i - 1, 0] = place_c(pts_b[i - 1, 0], pts_b[i - 1, 1])
        pts_c[i - 1, 1] = Decimal(0)
        pts_b[i, 0], pts_b[i, 1] = place_b(pts_a[i, 0], pts_a[i, 1], pts_c[i - 1, 0])
    return pts_b, pts_c


def find_periodic_branch(
        place_a: Callable[[int, Decimal, Core], np.ndarray], k_layers: int, bisection_precision: Decimal,
        core: Core) -> (np.ndarray, np.ndarray, np.ndarray, Decimal):
    r"""
    Carry out a bisection search of the parameter delta that influences the placement of the A disks on the convex
    polygonal chain so that a Böröczky packing with periodic boundary conditions is achieved.

    We parameterize the convex polygonal chains by the parameter delta that gives the distance of a line g_2^< below
    g_2 = 2 \sqrt{3}. All A disks are guaranteed to lie between g_3 = 2 + \sqrt{3} and g_2^<.

    For a given delta, this function constructs the A, B and C disks in the upper half plane in the (half-)branch up
    to the boundary of the system box for which its symmetry axis along the branch coincides with the positive x-axis.
    Here, it constructs k_layers A disks, k_layers B disks, and k_layers - 1 C disks.

    In order to allow for periodic boundary conditions, the x-distance between the disks A_{k_layers} and B_{k_layers}
    has to be 1 so that the branch has an additional symmetry axis at x = L / 2 going through B_{k_layers}. This only
    succeeds, if the A disks were placed on a convex polygonal chain with a correct curvature.

    This function carries out a bisection search of delta until the x-distance between A_{k_layers} and B_{k_layers}
    lies in the interval [1, 1 + bisection_precision]. It then returns the positions of the A, B, and C disks together
    with the final value of delta.

    Parameters
    ----------
    place_a : Callable[[int, Decimal, Core], np.ndarray[dtype=Decimal]]
        A callable that places the A disks on a convex polygonal chain between g_3 and g_2^< based on the number of
        layers, delta, and the core.
    k_layers : int
        The number of layers in the branch.
    bisection_precision : Decimal
        The precision of the bisection search.
    core : Core
        The core of the Böröczky packing.

    Returns
    -------
    (np.ndarray[dtype=Decimal], np.ndarray[dtype=Decimal], np.ndarray[dtype=Decimal], Decimal)
        The positions of the A disks, the positions of the B disks, the position of the C disks,
        the final value of delta.
    """
    min_delta, max_delta, distance_ab_x = Decimal(0), Decimal(1), Decimal(0)
    original_min_delta, original_max_delta = min_delta, max_delta
    while True:
        # For the A_circle, we actually vary the radius.
        delta = (min_delta + max_delta) / Decimal(2)
        pts_a = place_a(k_layers, delta, core)
        pts_b, pts_c = place_bc(k_layers, pts_a)
        distance_ab_x = pts_b[-1, 0] - pts_a[-1, 0]
        if distance_ab_x < Decimal(1):
            min_delta = delta
        elif distance_ab_x > Decimal(1) + bisection_precision:
            max_delta = delta
        else:
            break
        assert max_delta - min_delta > bisection_precision
    if min_delta == original_min_delta or max_delta == original_max_delta:
        warnings.warn("[WARNING] The true value of the distance delta that g_2^< lies below g_2 was probably not found "
                      "because one of the boundaries of the bisection search has still the initial value.",
                      RuntimeWarning)
    return pts_a, pts_b, pts_c, delta


def reflect_a_disks(pts_a: np.ndarray, core: Core) -> np.ndarray:
    """
    Reflect the A disks in one polygonal convex chain about all symmetry axes of the square to construct all A disks in
    the Böröczky packing.

    The given A disks should lie in the upper half plane in the (half-)branch up to the boundary of the system box whose
    symmetry axis along the branch coincides with the positive x-axis.

    The number of given A disks is equal to the number of layers k_layers of the branch. The number of all A disks in
    Böröczky packing depends on the core, because the A_1 disk is shared in different branches only for the Kahle core.

    For the Böröczky core, there are 8 * k_layers A disks after the reflections so that the shape of the returned numpy
    array is (8 * k_layers, 2).

    For the Kahle core, there are 8 * (k_layers - 1) + 4 A disks after the reflections so that the shape of the returned
    numpy array is (8 * (k_layers - 1) + 4, 2).

    Parameters
    ----------
    pts_a : np.ndarray[dtype=Decimal]
        The positions of the A disks in the upper half plan in the (half-)branch up to the boundary of the system box
        whose symmetry axis along the branch coincides with the positive x-axis.
    core : Core
        The core of the Böröczky packing.

    Returns
    -------
    np.ndarray[dtype=Decimal]
        The positions of all A disks in the Böröczky packing.
    """
    all_pts_a = []
    k_layers = pts_a.shape[0]
    if core is Core.BORO:
        lower_index = 0
    else:
        assert core is Core.KAHLE
        all_pts_a.append([pts_a[0, 0], pts_a[0, 1]])
        all_pts_a.append([-pts_a[0, 0], pts_a[0, 1]])
        all_pts_a.append([pts_a[0, 1], -pts_a[0, 0]])
        all_pts_a.append([-pts_a[0, 1], -pts_a[0, 0]])
        lower_index = 1
    for i in range(lower_index, k_layers):
        all_pts_a.append([-pts_a[i, 0], pts_a[i, 1]])
        all_pts_a.append([-pts_a[i, 0], -pts_a[i, 1]])
        all_pts_a.append([pts_a[i, 0], -pts_a[i, 1]])
        all_pts_a.append([pts_a[i, 0], pts_a[i, 1]])
        all_pts_a.append([pts_a[i, 1], -pts_a[i, 0]])
        all_pts_a.append([-pts_a[i, 1], pts_a[i, 0]])
        all_pts_a.append([-pts_a[i, 1], -pts_a[i, 0]])
        all_pts_a.append([pts_a[i, 1], pts_a[i, 0]])
    all_pts_a = np.array(all_pts_a)
    assert all_pts_a.shape == ((8 * k_layers, 2) if core is Core.BORO else (8 * (k_layers - 1) + 4, 2))
    return all_pts_a


def reflect_b_disks(pts_b: np.ndarray) -> np.ndarray:
    """
    Reflect the given B disks about all symmetry axes of the square to construct all A disks in the Böröczky packing.

    The given B disks should lie in the upper half plane in the (half-)branch up to the boundary of the system box whose
    symmetry axis along the branch coincides with the positive x-axis.

    The number of given B disks is equal to the number of layers k_layers of the branch. Since the last B_{k_layers}
    disk lies on the boundary of the system box, there are 8 * (k_layers - 1) + 4 B disks after the reflections so that
    the shape of the returned numpy array is (8 * (k_layers - 1) + 4, 2).

    Parameters
    ----------
    pts_b : np.ndarray[dtype=Decimal]
        The positions of the B disks in the upper half plan in the (half-)branch up to the boundary of the system box
        whose symmetry axis along the branch coincides with the positive x-axis.

    Returns
    -------
    np.ndarray[dtype=Decimal]
        The positions of all B disks in the Böröczky packing.
    """
    all_pts_b = []
    k_layers = pts_b.shape[0]
    for i in range(k_layers - 1):
        all_pts_b.append([-pts_b[i, 0], pts_b[i, 1]])
        all_pts_b.append([-pts_b[i, 0], -pts_b[i, 1]])
        all_pts_b.append([pts_b[i, 0], -pts_b[i, 1]])
        all_pts_b.append([pts_b[i, 0], pts_b[i, 1]])
        all_pts_b.append([pts_b[i, 1], -pts_b[i, 0]])
        all_pts_b.append([-pts_b[i, 1], pts_b[i, 0]])
        all_pts_b.append([-pts_b[i, 1], -pts_b[i, 0]])
        all_pts_b.append([pts_b[i, 1], pts_b[i, 0]])
    all_pts_b.append([pts_b[-1, 0], pts_b[-1, 1]])
    all_pts_b.append([pts_b[-1, 0], -pts_b[-1, 1]])
    all_pts_b.append([pts_b[-1, 1], pts_b[-1, 0]])
    all_pts_b.append([-pts_b[-1, 1], pts_b[-1, 0]])
    all_pts_b = np.array(all_pts_b)
    assert all_pts_b.shape == (8 * (k_layers - 1) + 4, 2)
    return all_pts_b


def reflect_c_disks(pts_c: np.ndarray) -> np.ndarray:
    """
    Reflect the given C disks about all symmetry axes of the square to construct all C disks in the Böröczky packing.

    The given C disks should lie in the (half-)branch up to the boundary of the system box whose symmetry axis along the
    branch coincides with the positive x-axis.

    The number of given C disks is equal to the number of layers k_layers of the branch minus one. Since the C disks
    lie on the x-axis, there are 4 * (k_layers - 1) C disks after the reflections so that the shape of the returned
    numpy array is (4 * (k_layers - 1), 2).

    Parameters
    ----------
    pts_c : np.ndarray[dtype=Decimal]
        The positions of the C disks in the (half-)branch up to the boundary of the system box whose symmetry axis along
        the branch coincides with the positive x-axis.

    Returns
    -------
    np.ndarray[dtype=Decimal]
        The positions of all C disks in the Böröczky packing.
    """
    all_pts_c = []
    k_layers_minus_one = pts_c.shape[0]
    for i in range(k_layers_minus_one):
        all_pts_c.append([pts_c[i, 0], pts_c[i, 1]])
        all_pts_c.append([-pts_c[i, 0], pts_c[i, 1]])
        all_pts_c.append([pts_c[i, 1], -pts_c[i, 0]])
        all_pts_c.append([pts_c[i, 1], pts_c[i, 0]])
    all_pts_c = np.array(all_pts_c)
    assert all_pts_c.shape == (4 * k_layers_minus_one, 2)
    return all_pts_c


def check_packing(all_pts: np.ndarray, system_length: Decimal, distance_precision: Decimal) -> int:
    """
    Check that the given Böröczky packing has no overlaps within the given distance precision, and that the number of
    contacts for each disk is at least three. Return the total number of contacts within the hard-disk packing.

    Two disks overlap if their distance (that is possibly corrected for periodic boundary conditions) is smaller than
    2 - distance_precision.

    Two disks are in contact if their distance lies within the interval
    [2 - distance_precision, 2 + distance_precision].

    Parameters
    ----------
    all_pts : np.ndarray[dtype=Decimal]
        The positions of the disks in the Böröczky packing.
    system_length : Decimal
        The side length of the periodic square box of the Böröczky packing.
    distance_precision : Decimal
        The distance precision.

    Returns
    -------
    int
        The number of contacts of the Böröczky packing.

    Raises
    ------
    RuntimeError
        If there are overlaps between disks in the Böröczky packing.
        If any disk in the Böröczky packing has less than three contacts.
    """
    number_contacts = 0
    number_overlaps = 0
    for i in range(all_pts.shape[0]):
        number_contacts_i = 0
        for j in range(all_pts.shape[0]):
            if i == j:
                continue
            distance_x = all_pts[i, 0] - all_pts[j, 0]
            if distance_x >= system_length / Decimal(2):
                distance_x -= system_length
            elif distance_x < -system_length / Decimal(2):
                distance_x += system_length
            distance_y = all_pts[i, 1] - all_pts[j, 1]
            if distance_y >= system_length / Decimal(2):
                distance_y -= system_length
            elif distance_y < -system_length / Decimal(2):
                distance_y += system_length
            distance = (distance_x * distance_x + distance_y * distance_y).sqrt()
            if distance < Decimal(2) - distance_precision:
                print(f"[ERROR] Disks with indices {i} and {j} overlap at distance {distance}.", file=sys.stderr)
                number_overlaps += 1
            if distance <= Decimal(2) + distance_precision:
                number_contacts_i += 1
                # Avoid double counting.
                if j < i:
                    number_contacts += 1
        if number_contacts_i < 3:
            raise RuntimeError(f"[ERROR] Disk with index {i} has not at least three contacts (number of "
                               f"contacts={number_contacts_i}).")
    if number_overlaps != 0:
        raise RuntimeError(f"[ERROR] Detected overlaps between disks in the Böröczky packing.")
    return number_contacts


def main() -> None:
    """Generate the Böröczky packing based on the command line arguments."""
    print_start_message()
    args = parse_options(sys.argv[1:])

    if args.core == "boro":
        print("Core: Böröczky")
        core = Core.BORO
    else:
        print("Core: Kahle")
        assert args.core == "kahle"
        core = Core.KAHLE
    if args.chain == "geometric":
        print(f"Convex polygonal chain of the A disks: Geometric with attenuation parameter phi={args.phi}")
        phi = Decimal(args.phi)
        place_a = lambda k_layers, d, c: place_a_geometric(k_layers, d, c, phi)
    else:
        print(f"Convex polygonal chain of the A disks: Circle")
        assert args.chain == "circle"
        place_a = place_a_circle
    layers = args.layers
    print(f"Number of layers: {layers}")
    if args.file == "":
        if core is Core.BORO:
            file = f"Boro_Layers{layers}.txt"
        else:
            assert core is Core.KAHLE
            file = f"Kahle_Layers{layers}.txt"
    else:
        file = args.file
    print(f"Output filename: {file}")
    precision = args.precision
    getcontext().prec = args.precision
    print(f"Precision in places: {precision}")
    if args.bisection <= 0:
        bisection_precision = Decimal(f"1e-{precision - 3}")
        distance_precision = Decimal(f"1e-{precision - 5}")
    else:
        bisection_precision = Decimal(f"1e-{args.bisection}")
        distance_precision = Decimal(f"1e-{args.bisection - 2}")
    print(f"Precision of the bisection: {bisection_precision}")
    plot = args.plot
    print(f"Plot the Böröczky packing: {plot}")

    print()

    print("Constructing the core.")
    pts_core = place_core(core)

    print("Starting bisection search for the parameter delta in the convex polygonal chain that allows for periodic "
          "boundary conditions.")
    pts_a, pts_b, pts_c, delta = find_periodic_branch(place_a, layers, bisection_precision, core)
    print(f"Finished bisection search with the final value {delta} for the parameter delta that gives the distance "
          f"of g_2^< below g_2 (the convex polygonal chain lies between g_2^< and g_2).")

    print("Reflecting the branch about the symmetry axes of the square box.")
    all_pts_a = reflect_a_disks(pts_a, core)
    all_pts_b = reflect_b_disks(pts_b)
    all_pts_c = reflect_c_disks(pts_c)

    all_pts = np.concatenate((all_pts_a, all_pts_b, all_pts_c, pts_core))
    system_length = Decimal(2) * pts_b.max()
    print(f"Finished construction of a periodic Böröczky packing in a square box with side length {system_length}.")
    # Correct the positions so that the square from 0 to L in each direction.
    all_pts += system_length / Decimal(2)

    if plot:
        print(f"Plotting Böröczky packing.")
        import matplotlib.pyplot as plt
        plt.figure()
        plt.xlim(0, system_length)
        plt.ylim(0, system_length)
        ax = plt.gca()
        for i in range(len(all_pts)):
            circle = plt.Circle(all_pts[i, :], 1)
            ax.add_patch(circle)
            ax.annotate(i, all_pts[i, :], va="center", ha="center")
        ax.set_aspect("equal")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
        plt.close()

    print(f"Checking the hard-disk packing with a precision of {distance_precision} (two orders of magnitudes larger "
          f"than the precision of the bisection search for delta).")
    number_of_contacts = check_packing(all_pts, system_length, distance_precision)
    if core is core.BORO:
        assert number_of_contacts == 32 * layers + 20
    else:
        assert core is core.KAHLE
        assert number_of_contacts == 32 * layers + 4
    print(f"Found {number_of_contacts} contacts in the Böröczky packing.")

    print(f"Writing Böröczky packing into the file {file}.")
    with open(file, "w") as f:
        print(f"# Core: {args.core}", file=f)
        print(f"# Chain: {args.chain}", file=f)
        print(f"# Phi: {args.phi}", file=f)
        print(f"# Layers: {args.layers}", file=f)
        print(f"# Precision: {args.precision}", file=f)
        print(f"# Bisection precision: {bisection_precision}", file=f)
        print(f"# Final value of delta after bisection search: {delta}", file=f)
        print(f"# Distance precision: {distance_precision}", file=f)
        print(f"# Number of contacts: {number_of_contacts}", file=f)
        print(f"# System length: {system_length}", file=f)
        for i in range(all_pts.shape[0]):
            print(f"{all_pts[i, 0]} {all_pts[i, 1]}", file=f)


if __name__ == '__main__':
    main()
