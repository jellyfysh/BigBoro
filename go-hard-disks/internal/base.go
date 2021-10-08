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
	"math/big"
)

// Return the squared euclidean norm of the given vector.
func normSquared(vector []*big.Float) *big.Float {
	sum := new(big.Float)
	mul := new(big.Float)
	for _, value := range vector {
		sum.Add(sum, mul.Mul(value, value))
	}
	return sum
}

// Return the dot product of the two given vectors.
func dot(vectorOne, vectorTwo []*big.Float) *big.Float {
	sum := new(big.Float)
	mul := new(big.Float)
	for i := range vectorOne {
		sum.Add(sum, mul.Mul(vectorOne[i], vectorTwo[i]))
	}
	return sum
}

// Return shortest separation vector with respect to boundary conditions in a cubic box with the given system length.
func separationVector(positionOne, positionTwo []*big.Float, sys *systemParameters) []*big.Float {
	separation := make([]*big.Float, len(positionOne))
	for i := range positionOne {
		separation[i] = new(big.Float)
		separation[i].Sub(positionTwo[i], positionOne[i])
		if separation[i].Cmp(sys.systemLengthOverTwo) >= 0 {
			separation[i].Sub(separation[i], sys.systemLength)
		} else if separation[i].Cmp(sys.negativeSystemLengthOverTwo) < 0 {
			separation[i].Add(separation[i], sys.systemLength)
		}
	}
	return separation
}

// Return shortest distance squared with respect to boundary conditions in a cubic box with the given system length.
func distance(positionOne, positionTwo []*big.Float, sys *systemParameters) *big.Float {
	distanceSquared := new(big.Float)
	separation := new(big.Float)
	for i := range positionOne {
		separation.Sub(positionTwo[i], positionOne[i])
		if separation.Cmp(sys.systemLengthOverTwo) >= 0 {
			separation.Sub(separation, sys.systemLength)
		} else if separation.Cmp(sys.negativeSystemLengthOverTwo) < 0 {
			separation.Add(separation, sys.systemLength)
		}
		distanceSquared.Add(distanceSquared, separation.Mul(separation, separation))
	}
	return distanceSquared.Sqrt(distanceSquared)
}
