# Event-chain Monte Carlo simulation of hard disks


The arbitrary-precision Go application `go-hard-disks` implements several variants of the event-chain Monte Carlo 
algorithm (ECMC) for two-dimensional hard-disk systems in a periodic square box. Straight, reflective, forward, and 
Newtonian ECMC are implemented (see, e.g., Section 3.1.2 of [[Hoellmer2021]](https://arxiv.org/abs/2109.13343) and 
references therein). 

## Installing

First, you need to install Go (see, e.g., https://golang.org/doc/install). We expect the Go application `go-hard-disks` 
to work with any Go version >= 1.13. We tested it with Go 1.16.

Afterwards, open a terminal and go into the [`cmd/go-hard-disks/`](cmd/go-hard-disks) directory. There you can use 
one of the following options for the installation of the application:
1. The command `go install` installs a `go-hard-disks` executable to the directory named by the `GOBIN` environment 
variable, which defaults to `$GOPATH/bin/` or `$HOME/go/bin/` if the `GOPATH` environment variable is not set. You might 
want to add the `GOBIN` directory to your `PATH` environment variable for convenience.
2. The command `go build` installs a `go-hard-disks` executable to the current directory (that is, the
[`cmd/go-hard-disks/`](cmd/go-hard-disks) directory).
3. The command `go build -o output` installs the executable with the name and path specified in `output` (that is, 
`go build -o ../../ecmc-hard-disks` creates a `ecmc-hard-disks` executable in the directory of this `README.md`).

## Using

All parameters of the hard-disk ECMC simulation are specified by command-line arguments. In order to get an overview
of the possible command-line arguments, use the `--help` command-line option. If the `go-hard-disks` executable is,
for example, installed in the `~/go/bin/` directory (which is the default path of `go install` after a new installation 
of Go), use
```
~/go/bin/go-hard-disks --help
```
The possible arguments and their default values are as follows:

1. `-ecmc string`: specify the ECMC scheme; possible values are 'straight', 'reflective', 'forward', and 'newtonian'
   (default "forward")
2. `-delphi float`: specify the rotation angle in degrees of the direction of straight ECMC after each chain time
(ignored if reflective, forward, or Newtonian ECMC scheme is specified); if 0.0, the direction alternates between the
positive x- and y-directions; if negative, the direction is sampled randomly after each chain time (default 0)
3. `-chain float`: specify the chain time of ECMC after which the active disk and its velocity are resampled; if 0.0 or
   negative, the chain time is set to infinity (default 0)
4. `-sampling-t float`: specify the sampling time after which the maximum nearest-neighbor distance is sampled 
(default 0.1)
5. `-sampling-n uint`: specify the total number of samples after which the simulation is finished (default 1000)
6. `-init string`: specify the filename of the initial configuration (default "Initial.txt")
7. `-final string`: specify the filename for the storage of the final configuration; if empty, the final configuration 
is not stored (default '')
8. `-file string`: specify the output filename that stores the sampled maximum nearest-neighbor distances 
(default "MaximumNearestNeighborDistance.dat")
9. `-length string`: specify the system length of the periodic square box (default "25.451518000277975")
10. `-sigma string`: specify the hard-disk radius sigma (default "1.0")
11. `-cells int`: specify the number cells per side in the cell-occupancy system (default 12)
12. `-n int`: specify the number of hard disks (default 96)
13. `-rattle string`: specify the radius of the circle in which the initial positions of the hard-disks are rattled; if 
0.0 or negative, the original initial positions are used (default "0.0")
14. `-prec uint`: specify the number of mantissa bits that are used for arbitrary-precision arithmetic (default 53)
15. `-max-events uint`: specify the maximum number of events after which the simulation is interrupted 
(default 18446744073709551615)

The command-line arguments 9., 10., and 13. have the type `string` in the above list. These strings are directly parsed 
into a multi-precison floating point number (with the number of mantissa bits specified by the `-prec` argument) which 
circumvents loosing precision by first parsing the value into a standard double. The given string values should be text 
representations of floats in base 10, and are passed to the `Parse` function of Go's `math/big` package 
(see https://pkg.go.dev/math/big#Float.Parse).

The initial hard-disk configuration file should contain the x- and y-positions of the hard disks, separated by a space 
or a tab, in each non-header row. Ignored header rows can be used for comments and start with `#`. Such a file format 
is, for example, the output of the Python script 
[`construct_packing.py`](../construct_packing/construct_packing.py) of the BigBoro software package. The `-n` 
command-line argument is solely used to assert that the right number of hard-disk positions was parsed from the file of 
the initial configuration. The final configuration is (optionally) stored in a file with the same format.

If the specified hard-disk radius in the `-sigma` command-line arguments yields overlapping hard disks in the initial 
configuration (under consideration of periodic boundary conditions in a square box with the specified system length), 
the application will panic. The changing configuration over the course of the simulation is likewise checked each 
10000 events, 1000 samples, and at the end of the simulation.

The application relies on a cell-occupancy system for the efficient simulation of systems with a large number of hard 
disks. The minimum number of cells per side is three. Each cell should be bigger than two times the hard-disk radius. 
This sets an upper bound on the number of cells per side. If the specified number of cells per side falls outside the 
allowed region, the application will panic.

The application samples the maximum nearest-neighbor distance of the hard disks after each sampling time (see eq. (8) of
[[Hoellmer2021]](https://arxiv.org/abs/2109.13343)). The specified sampling file stores the sampled maximum 
nearest-neighbor distance together with the running number of events in each row.

## Examples

As an example, we consider the Böröczky packing with the Kahle core and a 5-layered geometric convex polygonal chain 
with attenuation parameter phi=0.7. This packing can be generated with the Python script 
[`construct_packing.py`](../construct_packing/construct_packing.py) of the BigBoro software package, and is also part of 
its example packings (see 
[`../construct_packing/example_packings/kahle_geometric_5/packing.txt`](../construct_packing/example_packings/kahle_geometric_5/packing.txt)).
Copy the example packing into a directory in which you want to run the ECMC simulation.

We consider the escape of forward ECMC from an epsilon-relaxed Böröczky configuration, where the original hard-disk 
radius is multiplied by (1 - epsilon) with epsilon=10^{-10}. We use an infinite chain time. As the precision, we choose 
200 mantissa bits. Under the assumption that you installed the Go application `go-hard-disks` using `go install` to 
create the executable `~/go/bin/go-hard-disks`, execute the following command to sample 100 maximum nearest-neighbor 
distances and the corresponding number of events in intervals of the sampling time 0.1 in the file 
`MaximumNearestNeighborDistances.txt`:
```
~/go/bin/go-hard-disks -cells=12 -chain=0.0 -ecmc=forward -file=MaximumNearestNeighborDistance.dat -init=packing.txt -n=96 -prec=200 -sampling-n=100 -sampling-t=0.1 -sigma=0.9999999999 -length=25.451518000277973442489180458890531734216606953961422451578118800844048025977126841626591893011618209821314411437263008910423767729009689705721865596107365845523336502309341340960867419256611371168060
```

## Remarks

- In the Go application `go-hard-disks`, the central simulation box ranges from 0 to the system length L in each 
direction (in contrast to [[Hoellmer2021]](https://arxiv.org/abs/2109.13343) where it ranges from -L/2 to L/2).
- In its current form, the application only samples the maximum nearest-neighbor distance and the number of events. 
It can be extended to other observables by modifying the `Write` function in line 57 of 
[`internal/file.go`](internal/file.go).
- Standard double or single float precision are achieved by setting the number of mantissa bits to 53 or 24.
- Although the application implements a cell-occupancy system, ECMC simulations of a large number of hard disks can 
become slow. This is because the frequent checks for overlaps between hard disks do not rely on the cell-occupancy 
system at the moment. Moreover, the maximum nearest-neighbor distance is also computed without the help of the 
cell-occupancy system. This might change in future versions.
