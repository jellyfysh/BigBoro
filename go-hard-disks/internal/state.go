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

package internal

import (
	"bufio"
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"os"
	"strings"
)

type Sphere struct {
	position, velocity []*big.Float
}

type Spheres []Sphere

type identifier int

type systemParameters struct {
	prec                                                                                              uint
	systemLength, systemLengthOverTwo, negativeSystemLengthOverTwo, sigma, twoSigma, fourSigmaSquared *big.Float
}

func SystemParameters(L, s string, prec uint) *systemParameters {
	systemLength, _, err := new(big.Float).SetPrec(prec).Parse(L, 10)
	if err != nil {
		panic(fmt.Sprintf("Expected a float as the system length: %s", L))
	}
	sigma, _, err := new(big.Float).SetPrec(prec).Parse(s, 10)
	if err != nil {
		panic(fmt.Sprintf("Expected a float as the radius sigma: %s", s))
	}

	two := new(big.Float).SetPrec(prec).SetInt64(2)
	systemLengthOverTwo := new(big.Float)
	negativeSystemLengthOverTwo := new(big.Float)
	twoSigma := new(big.Float)
	fourSigmaSquared := new(big.Float)

	systemLengthOverTwo.Quo(systemLength, two)
	negativeSystemLengthOverTwo.Neg(systemLengthOverTwo)
	twoSigma.Mul(two, sigma)
	fourSigmaSquared.Mul(twoSigma, twoSigma)
	sys := systemParameters{
		prec,
		systemLength,
		systemLengthOverTwo,
		negativeSystemLengthOverTwo,
		sigma,
		twoSigma,
		fourSigmaSquared}
	return &sys
}

func computeMaximumNearestNeighborDistance(spheres Spheres, sys *systemParameters) *big.Float {
	var dist *big.Float
	var sphereOne *Sphere
	var sphereIndexTwo int
	nearestNeighborDistance := new(big.Float)
	maxNearestNeighborDistance := new(big.Float).SetPrec(sys.prec)
	tenSystemLengths := new(big.Float).SetPrec(sys.prec).SetInt64(10)
	tenSystemLengths.Mul(tenSystemLengths, sys.systemLength)
	for sphereIndexOne := range spheres {
		sphereOne = &spheres[sphereIndexOne]
		nearestNeighborDistance.Copy(tenSystemLengths)
		for sphereIndexTwo = range spheres {
			if sphereIndexOne != sphereIndexTwo {
				dist = distance(sphereOne.position, spheres[sphereIndexTwo].position, sys)
				if dist.Cmp(nearestNeighborDistance) < 0 {
					nearestNeighborDistance.Copy(dist)
				}
			}
		}
		if nearestNeighborDistance.Cmp(tenSystemLengths) == 0 {
			panic("Error while computing maximum nearest-neighbor distance.")
		}
		if nearestNeighborDistance.Cmp(maxNearestNeighborDistance) > 0 {
			maxNearestNeighborDistance.Copy(nearestNeighborDistance)
		}
	}
	return maxNearestNeighborDistance
}

func ReadSpheresTxt(numberSpheres int, sys *systemParameters, filename string, rattleRadiusString string) Spheres {
	spheres := make(Spheres, numberSpheres)
	for i := range spheres {
		spheres[i].position = make([]*big.Float, 2)
		spheres[i].velocity = nil
	}

	// Open input file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	// close input file on exit and check for its returned error
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()
	scanner := bufio.NewScanner(file)
	scanner.Split(bufio.ScanLines)
	var line string
	var splitLine []string
	currentSphereIndex := 0

	rattleRadius, _, err := new(big.Float).SetPrec(sys.prec).Parse(rattleRadiusString, 10)
	if err != nil {
		panic(fmt.Sprintf("Expected a float as the rattle radius: %s", rattleRadiusString))
	}

	for scanner.Scan() {
		line = scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		if currentSphereIndex >= numberSpheres {
			panic("Too many spheres in file.")
		}
		splitLine = strings.Split(line, " ")
		if len(splitLine) != 2 {
			splitLine = strings.Split(line, "\t")
			if len(splitLine) != 2 {
				panic(fmt.Sprintf("Unexpected line: %s", line))
			}
		}
		spheres[currentSphereIndex].position[0], _, err = new(big.Float).SetPrec(sys.prec).Parse(splitLine[0], 10)
		if err != nil {
			panic(fmt.Sprintf("Expected float in first column of line: %s", line))
		}
		spheres[currentSphereIndex].position[1], _, err = new(big.Float).SetPrec(sys.prec).Parse(splitLine[1], 10)
		if err != nil {
			panic(fmt.Sprintf("Expected float in second column of line: %s", line))
		}

		if rattleRadius.Sign() > 0 {
			randomAngle := -math.Pi + 2.0*math.Pi*rand.Float64()
			cosAngle := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomAngle))
			sinAngle := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomAngle))
			randFloat := new(big.Float).SetPrec(sys.prec).SetFloat64(rand.Float64())
			randRadius := new(big.Float).SetPrec(sys.prec)
			// Generate random point in circle. See: https://stackoverflow.com/questions/5837572
			randRadius.Mul(rattleRadius, randFloat.Sqrt(randFloat))
			t := new(big.Float)
			spheres[currentSphereIndex].position[0].Add(spheres[currentSphereIndex].position[0],
				t.Mul(randRadius, cosAngle))
			spheres[currentSphereIndex].position[1].Add(spheres[currentSphereIndex].position[1],
				t.Mul(randRadius, sinAngle))
		}

		if spheres[currentSphereIndex].position[0].Cmp(sys.systemLength) >= 0 {
			spheres[currentSphereIndex].position[0].Sub(spheres[currentSphereIndex].position[0], sys.systemLength)
		} else if spheres[currentSphereIndex].position[0].Sign() < 0 {
			spheres[currentSphereIndex].position[0].Add(spheres[currentSphereIndex].position[0], sys.systemLength)
		}
		if spheres[currentSphereIndex].position[1].Cmp(sys.systemLength) >= 0 {
			spheres[currentSphereIndex].position[1].Sub(spheres[currentSphereIndex].position[1], sys.systemLength)
		} else if spheres[currentSphereIndex].position[1].Sign() < 0 {
			spheres[currentSphereIndex].position[1].Add(spheres[currentSphereIndex].position[1], sys.systemLength)
		}
		currentSphereIndex += 1
	}
	if currentSphereIndex != numberSpheres {
		panic(fmt.Sprintf("The file %s contains %d spheres, but %d spheres were specified in the command line "+
			"argument.", filename, currentSphereIndex, numberSpheres))
	}
	checkState(spheres, sys, true)
	return spheres
}

func checkState(spheres Spheres, sys *systemParameters, exact bool) {
	var dist *big.Float
	var sphereIndexTwo int
	var sphereOne, sphereTwo *Sphere
	bound := new(big.Float).SetPrec(sys.prec)
	if exact {
		bound.Copy(sys.twoSigma)
	} else {
		power := int(float64(sys.prec) * math.Log10(2.0))
		small := new(big.Float).SetPrec(sys.prec)
		_, _, err := small.Parse(fmt.Sprintf("1.0e-%d", power-2), 10)
		if err != nil {
			panic(fmt.Sprintf("Could not parse precision string for check: 1.0e-%d", power-2))
		}
		bound.Sub(sys.twoSigma, small)
	}
	for sphereIndexOne := range spheres {
		sphereOne = &spheres[sphereIndexOne]
		for sphereIndexTwo = sphereIndexOne + 1; sphereIndexTwo < len(spheres); sphereIndexTwo++ {
			sphereTwo = &spheres[sphereIndexTwo]
			dist = distance(sphereOne.position, sphereTwo.position, sys)
			if dist.Cmp(bound) <= 0 {
				fmt.Println("Problem between sphere", sphereIndexOne, "and", sphereIndexTwo, ".")
				fmt.Println("Position sphere one:", sphereOne.position)
				fmt.Println("Position sphere two:", sphereTwo.position)
				fmt.Printf("Distance: %v\n", dist)
				panic("Illegal state detected.")
			}
		}
	}
}

func WriteStateToFile(spheres Spheres, sys *systemParameters, filename string) {
	// Open output file
	file, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	// close output file on exit and check for its returned error
	defer func() {
		if err := file.Close(); err != nil {
			panic(err)
		}
	}()

	if _, err = file.WriteString(fmt.Sprintf("# N=%d, L=%v, sigma=%v\n",
		len(spheres), sys.systemLength, sys.sigma)); err != nil {
		panic(err)
	}

	if _, err = file.WriteString("# x y\n"); err != nil {
		panic(err)
	}

	power := int(float64(sys.prec) * math.Log10(2.0))
	for sphereIndex := range spheres {
		if _, err = file.WriteString(fmt.Sprintf("%s %s\n",
			spheres[sphereIndex].position[0].Text('f', power),
			spheres[sphereIndex].position[1].Text('f', power))); err != nil {
			panic(err)
		}
	}
}
