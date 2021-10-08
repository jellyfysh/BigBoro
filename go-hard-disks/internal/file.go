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
	"fmt"
	"os"
)

type File interface {
	Write(spheres Spheres, sys *systemParameters, numberEvents uint64)

	Finished() bool

	Close()
}

type MaximumNearestNeighborDistanceFile struct {
	osFile        *os.File
	numberSamples uint64
	takenSamples  *uint64
}

func SetUpMaximumNearestNeighborDistanceFile(filename string, numberSamples uint64) *MaximumNearestNeighborDistanceFile {
	osFile, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	if _, err = osFile.WriteString("# Maximum nearest-neighbor distance\n"); err != nil {
		panic(err)
	}
	fmt.Println("Writing maximum nearest-neighbor distance into the file", filename+".")
	var takenSamples uint64 = 0
	file := MaximumNearestNeighborDistanceFile{osFile, numberSamples, &takenSamples}
	return &file
}

func (file *MaximumNearestNeighborDistanceFile) Write(spheres Spheres, sys *systemParameters,
	numberEvents uint64) {

	maximumNearestNeighborDistance := computeMaximumNearestNeighborDistance(spheres, sys)
	if _, err := file.osFile.WriteString(
		fmt.Sprintf("%v\t%d\n", maximumNearestNeighborDistance, numberEvents)); err != nil {
		panic(err)
	}

	*file.takenSamples += 1
	if *file.takenSamples%1000 == 0 {
		checkState(spheres, sys, false)
		fmt.Println("Taken", *file.takenSamples, "samples.")
	}
}

func (file *MaximumNearestNeighborDistanceFile) Finished() bool {
	return *file.takenSamples >= file.numberSamples
}

func (file *MaximumNearestNeighborDistanceFile) Close() {
	if err := file.osFile.Close(); err != nil {
		panic(err)
	}
}
