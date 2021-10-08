# Collective escape modes


The double-precision Python script [`collective_escape_modes.py`](collective_escape_modes.py) identifies the orthonormal 
basis vectors Delta_a of the escape matrix M^{esc} from a two-dimensional packing 
x = {x_1, y_1, x_2, y_2, ..., x_N, y_N} of N hard disks in a periodic square box that have zero singular values (see 
Section 2.2.2 of [[Hoellmer2021]](https://arxiv.org/abs/2109.13343)). This is the solution space for 2N-dimensional 
displacements Delta = {Delta^x_1, Delta^y_1, \Delta^x_2, \Delta^y_2, ...}, where Delta^x_1 is the displacement of hard 
disk 1 along x, etc. Displacements from this solution space leave all contacts intact up to first order in the 
displacements Delta_i, and increases separations between disks in contact up to second order (see eq. (5)). They are 
thus collective infinitesimal displacements of all disks that escape from the packing.

## Installing

The Python script [`collective_escape_modes.py`](collective_escape_modes.py) can be executed with any Python3 
implementation. Its only required external dependency is [NumPy](https://numpy.org). An optional plotting part 
requires [Matplotlib](https://matplotlib.org). Both packages can be installed, e.g., using `pip` or `conda` (see
https://numpy.org/install/ and https://matplotlib.org/stable/users/installing.html).

We expect the Python script [`collective_escape_modes.py`](collective_escape_modes.py) to work with any Python3 version 
and any NumPy version. We tested it with cPython 3.9 and NumPy 1.21.

## Using

All parameters for the computation of the collective escape modes are specified by command-line arguments. For an overview over the available command-line arguments, use:
```
python3 collective_escape_modes.py --help
```

The filename for storing the hard-disk packing is the only required positional argument.

Further (optional) command-line arguments and their default values are as follows:

1. `--length LENGTH`: specify the system length of the periodic square box; if None, try to parse the system length from 
the file that stores the packing (default=None)
2. `--file FILE`: set the output filename for the collective escape modes; if empty, the escape modes are stored in 
'EscapeModes.txt' (default ='')
3. `--plot PLOT`: set the pdf-filename for the plots of all escape modes (requires matplotlib); if empty, no plot is 
created (default='')
4. `--version`: show program's version number and exit
5. `-h, --help`: show the help message and exit

The hard-disk packing file should contain the x- and y-positions of the hard disks, separated by any whitespace, in each 
non-header row. Header rows can be used for comments and start with `#`. Such a file format is, for example, the output 
of the Python script [`construct_packing.py`](../construct_packing/construct_packing.py) of the BigBoro software 
package.

If the `--length` command-line argument is omitted, the Python script 
[`collective_escape_modes.py`](collective_escape_modes.py) expects in the packing file a header line of the format 
'# System length: {system_length}' where {system_length} is a float. Such a line is, for example, included by the Python 
script [`construct_packing.py`](../construct_packing/construct_packing.py) of the BigBoro software package.

The script first considers all pairs of hard disks at distances in the interval [2 - 10^{-13}, 2 + 10^{-13}] as 
contacts (and errors if no contacts are found). Following the singular value decomposition, the script asserts that the
configurations x + 10^{-8} Delta_a for every orthonormal basis vector Delta_a are without overlaps and that all contacts 
persist at a precision 10^{-8}. Furthermore, the script ensures that the uniform translations of all disks along the x- 
and y-axis are part of the solution space. 

The basis vectors Delta_a are stored in each (non-header) line in a human-readable format in the specified output file 
of the `--file` command-line argument. The optional plotting part represents the basis vectors Delta_a as arrows 
starting from the centers of the hard disks (compare Fig. 2 of [[Hoellmer2021]](https://arxiv.org/abs/2109.13343)).

## Examples

Executing the following command generates the file `CollectiveEscapeModes.pdf` (with the help of Matplotlib):
```
python3 collective_escape_modes.py --plot=CollectiveEscapeModes.pdf ../construct_packing/example_packings/kahle_geometric_5/packing.txt 
```
The escape modes 25 to 27 (on the pages 26 to 28) are displayed in Fig. 2 of 
[[Hoellmer2021]](https://arxiv.org/abs/2109.13343). (Fig. 2 was generated with Python 3.9.6 and NumPy 1.21.2. It is 
possible that different versions might yield different basis vectors Delta'_a. It is only guarantueed that all basis vectors 
Delta_a belong to the solution space spanned by Delta'_a.)

## Remarks

- The final checks of the configurations x + 10^{-8} Delta_a for possible (wrong) overlaps and too few contacts 
currently compute the distance between every pair of hard disks. These checks thus becomes very slow for large packings. 
- Because of floating point arithmetic, it is possible that some singular values that vanish in theory are reported as 
nonzero. The script currently assumes that all singular values smaller than 10^{-13} are effectively zero, and prints a 
warning if any nonzero singular values are set to zero. The user should check carefully that the corresponding basis 
vectors are indeed part of the solution space for collective escape modes.
