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
	"math/big"
)

func SimulateWithCells(spheres Spheres, systemParameters *systemParameters, samplingTime float64,
	chainTime float64, cells *Cells, velocityUpdater VelocityUpdater, file File, maximumNumberEvents uint64) {

	var numberEvents uint64 = 0
	actId := velocityUpdater.initializeVelocities(spheres, systemParameters)
	chainTimeBig := new(big.Float).SetPrec(systemParameters.prec)
	if chainTime <= 0.0 {
		chainTimeBig.SetInf(false)
	} else {
		chainTimeBig.SetFloat64(chainTime)
	}
	samplingTimeBig := new(big.Float).SetPrec(systemParameters.prec).SetFloat64(samplingTime)
	remainingSamplingTime := new(big.Float).SetPrec(systemParameters.prec).Copy(samplingTimeBig)

	for !file.Finished() {
		actId = runChainCells(spheres, actId, chainTimeBig, systemParameters, cells, velocityUpdater, &numberEvents,
			samplingTimeBig, remainingSamplingTime, file, maximumNumberEvents)
		if actId == -1 {
			fmt.Println("Interrupted simulation with", numberEvents, "events.")
			break
		}
		actId = velocityUpdater.updateVelocitiesEndOfChain(spheres, actId, systemParameters)
	}

	checkState(spheres, systemParameters, false)
}
