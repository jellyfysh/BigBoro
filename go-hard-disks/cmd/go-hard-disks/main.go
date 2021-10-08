// BigBoro - Arbitrary-precision Python/Go software for Böröczky packings with ECMC computations
// https://github.com/jellyfysh/BigBoro
// Copyright (C) 2021 Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, and Werner Krauth
//
// This file is part of BigBoro.
//
// BigBoro is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// BigBoro is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with BigBoro in the LICENSE file.
// If not, see <https://www.gnu.org/licenses/>.
//
// If you use BigBoro in published work, please cite the following reference (see [Hoellmer2021] in References.bib):
// Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, and Werner Krauth,
// Sparse hard-disk packings and Markov chains,
// arXiv e-prints: 2109.13343 (2021), https://arxiv.org/abs/2109.13343.
//

package main

import (
	"flag"
	"fmt"
	"github.com/jellyfysh/BigBoro/go-hard-disks/internal"
	"math"
	"math/rand"
	"time"
)

func main() {
	fmt.Println("BigBoro (version 1.0) - Arbitrary-precision Python/Go software for Böröczky packings with " +
		"ECMC computations - https://github.com/jellyfysh/BigBoro")
	fmt.Println("Copyright (C) 2021 Philipp Höllmer, Nicolas Noirault, Botao Li, A. C. Maggs, and Werner Krauth")
	fmt.Println("")

	seed := time.Now().UnixNano()
	defer func() {
		if err := recover(); err != nil {
			fmt.Printf("Program started with random seed %d panicked.\n", seed)
			panic(err)
		}
	}()
	rand.Seed(seed)

	var chainTime, deltaPhiDegree, samplingTime float64
	var cellsPerSide, N int
	var numberSamples, maximumNumberEvents uint64
	var filename, initialConfiguration, finalConfiguration, ecmcScheme, systemLength, sigma, rattleRadius string
	var prec uint

	flag.Float64Var(&chainTime, "chain", 0.0,
		"specify the chain time of ECMC after which the active disk and its velocity are resampled; " +
		"if 0.0 or negative, the chain time is set to infinity (default 0)")
	flag.StringVar(&systemLength, "length", "25.451518000277975",
		"specify the system length of the periodic square box")
	flag.StringVar(&sigma, "sigma", "1.0", "specify the hard-disk radius sigma")
	flag.Float64Var(&deltaPhiDegree, "delphi", 0.0,
		"specify the rotation angle in degrees of the direction of straight ECMC after each chain time " +
		"(ignored if reflective, forward, or Newtonian ECMC scheme is specified); " +
		"if 0.0, the direction alternates between the positive x- and y-directions; " +
		"if negative, the direction is sampled randomly after each chain time (default 0)")
	flag.Float64Var(&samplingTime, "sampling-t", 0.1,
		"specify the sampling time after which the maximum nearest-neighbor distance is sampled")
	flag.Uint64Var(&numberSamples, "sampling-n", 1000,
		"specify the total number of samples after which the simulation is finished")
	flag.IntVar(&cellsPerSide, "cells", 12,
		"specify the number cells per side in the cell-occupancy system")
	flag.IntVar(&N, "n", 96,
		"specify the number of hard disks")
	flag.StringVar(&filename, "file", "MaximumNearestNeighborDistance.dat",
		"specify the output filename that stores the sampled maximum nearest-neighbor distances")
	flag.StringVar(&initialConfiguration, "init", "Initial.txt",
		"specify the filename of the initial configuration")
	flag.StringVar(&finalConfiguration, "final", "",
		"specify the filename for the storage of the final configuration; if empty, the final configuration is" +
		" not stored (default '')")
	flag.StringVar(&ecmcScheme, "ecmc", "forward",
		"specify the ECMC scheme; possible values are 'straight', 'reflective', 'forward', and 'newtonian'")
	flag.UintVar(&prec, "prec", 53,
		"specify the number of mantissa bits that are used for arbitrary-precision arithmetic")
	flag.Uint64Var(&maximumNumberEvents, "max-events", math.MaxUint64,
		"specify the maximum number of events after which the simulation is interrupted")
	flag.StringVar(&rattleRadius, "rattle", "0.0",
		"specify the radius of the circle in which the initial positions of the hard-disks are rattled; " +
		"if 0.0 or negative, the original initial positions are used")
	flag.Parse()
	fmt.Println("Ecmc scheme:", ecmcScheme)
	fmt.Println("Delta phi degree:", deltaPhiDegree)
	fmt.Println("Chain time:", chainTime)
	fmt.Println("Sampling time:", samplingTime)
	fmt.Println("Number of samples:", numberSamples)
	fmt.Println("Filename of initial configuration:", initialConfiguration)
	fmt.Println("Filename for final configuration:", finalConfiguration)
	fmt.Println("Filename for sampled values:", filename)
	fmt.Println("System length:", systemLength)
	fmt.Println("Sigma:", sigma)
	fmt.Println("Cells per side:", cellsPerSide)
	fmt.Println("Number of hard disks:", N)
	fmt.Println("Rattle radius for initial state:", rattleRadius)
	fmt.Println("Precision:", prec)
	fmt.Println("Maximum number of events:", maximumNumberEvents)

	var velocityUpdater internal.VelocityUpdater
	if ecmcScheme == "straight" {
		fmt.Println("Running straight ECMC.")
		if deltaPhiDegree == 0.0 {
			fmt.Println("Running periodic xy-scheme (since delta phi degree = 0.0).")
		} else if deltaPhiDegree < 0.0 {
			fmt.Println("Choosing direction randomly (since delta phi degree < 0.0).")
		} else {
			fmt.Println("Choosing direction sequentially with delta phi =", deltaPhiDegree, "deg.")
		}
		velocityUpdater = &internal.StraightVelocityUpdater{DeltaPhi: deltaPhiDegree * math.Pi / 180.0}
	} else if ecmcScheme == "reflective" {
		fmt.Println("Running reflective ECMC.")
		velocityUpdater = &internal.ReflectiveVelocityUpdater{}
	} else if ecmcScheme == "forward" {
		fmt.Println("Running forward ECMC.")
		velocityUpdater = &internal.ForwardVelocityUpdater{}
	} else if ecmcScheme == "newtonian" {
		fmt.Println("Running Newtonian ECMC.")
		beta := 1.0
		mass := 1.0
		velocityUpdater = &internal.NewtonianVelocityUpdater{StandardDeviation: math.Sqrt(1.0 / (beta * mass))}
	} else {
		panic("Ecmc scheme chosen by the -ecmc command line option can only be one of 'straight', 'reflective', " +
			"'forward', or 'newtonian'.")
	}

	systemParameters := internal.SystemParameters(systemLength, sigma, prec)
	fmt.Println("Reading initial state from .txt file.")
	spheres := internal.ReadSpheresTxt(N, systemParameters, initialConfiguration, rattleRadius)
	cells := internal.CreateCells(spheres, cellsPerSide, systemParameters)

	file := internal.SetUpMaximumNearestNeighborDistanceFile(filename, numberSamples)
	defer file.Close()

	internal.SimulateWithCells(spheres, systemParameters, samplingTime, chainTime, cells, velocityUpdater, file,
		maximumNumberEvents)

	if finalConfiguration != "" {
		fmt.Println("Writing final configuration.")
		internal.WriteStateToFile(spheres, systemParameters, finalConfiguration)
	}
}
