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

type event struct {
	time   *big.Float
	nextId identifier
}

func timeToSphereEvent(activeSphere, targetSphere *Sphere, sys *systemParameters) *big.Float {
	separation := separationVector(activeSphere.position, targetSphere.position, sys)
	separationSquared := normSquared(separation)
	velocityDotSeparation := dot(activeSphere.velocity, separation)
	velocitySquared := normSquared(activeSphere.velocity)
	if velocityDotSeparation.Sign() > 0 {
		sqrtTerm := new(big.Float)
		t := new(big.Float)
		t2 := new(big.Float)
		t.Sub(separationSquared, sys.fourSigmaSquared)
		t.Mul(t, velocitySquared)
		t2.Mul(velocityDotSeparation, velocityDotSeparation)
		sqrtTerm.Sub(t2, t)
		if sqrtTerm.Sign() > 0 {
			sqrtTerm.Sqrt(sqrtTerm)
			t.Sub(velocityDotSeparation, sqrtTerm)
			return t.Quo(t, velocitySquared)
		}
	}
	// Positive infinity.
	return new(big.Float).SetPrec(sys.prec).SetInf(false)
}
