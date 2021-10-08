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
	"math"
	"math/big"
	"math/rand"
)

type VelocityUpdater interface {
	initializeVelocities(spheres Spheres, sys *systemParameters) identifier

	updateVelocityEvent(activeSphere, targetSphere *Sphere, sys *systemParameters)

	updateVelocitiesEndOfChain(spheres Spheres, actId identifier, sys *systemParameters) identifier

	correctVelocitiesSampling(spheres Spheres, actId identifier, sys *systemParameters)
}

type StraightVelocityUpdater struct {
	DeltaPhi float64
	// See: https://stackoverflow.com/questions/44370277
	currentPhi *float64
}

func (u *StraightVelocityUpdater) initializeVelocities(spheres Spheres, sys *systemParameters) identifier {
	actId := identifier(rand.Intn(len(spheres)))
	if u.DeltaPhi == 0.0 {
		velocityX := new(big.Float).SetPrec(sys.prec).SetInt64(1)
		velocityY := new(big.Float).SetPrec(sys.prec)
		spheres[actId].velocity = []*big.Float{velocityX, velocityY}
	} else {
		randomPhi := -math.Pi + 2.0*math.Pi*rand.Float64()
		u.currentPhi = &randomPhi
		velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomPhi))
		velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomPhi))
		spheres[actId].velocity = []*big.Float{velocityX, velocityY}
	}
	return actId
}

func (u *StraightVelocityUpdater) updateVelocityEvent(activeSphere, targetSphere *Sphere, _ *systemParameters) {
	targetSphere.velocity = activeSphere.velocity
	activeSphere.velocity = nil
}

func (u *StraightVelocityUpdater) updateVelocitiesEndOfChain(spheres Spheres, actId identifier,
	sys *systemParameters) identifier {

	newActId := identifier(rand.Intn(len(spheres)))
	oldActive := &spheres[actId]
	newActive := &spheres[newActId]
	if newActId != actId {
		newActive.velocity = make([]*big.Float, 2)
	}
	if u.DeltaPhi == 0.0 {
		newActive.velocity[0], newActive.velocity[1] = oldActive.velocity[1], oldActive.velocity[0]
	} else if u.DeltaPhi < 0.0 {
		randomPhi := -math.Pi + 2.0*math.Pi*rand.Float64()
		u.currentPhi = &randomPhi
		velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomPhi))
		velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomPhi))
		newActive.velocity[0] = velocityX
		newActive.velocity[1] = velocityY
	} else {
		nextPhi := *u.currentPhi + u.DeltaPhi
		if nextPhi >= math.Pi {
			nextPhi -= 2.0 * math.Pi
		}
		u.currentPhi = &nextPhi
		velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(nextPhi))
		velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(nextPhi))
		newActive.velocity[0] = velocityX
		newActive.velocity[1] = velocityY
	}
	if newActId != actId {
		oldActive.velocity = nil
	}
	return newActId
}

func (u *StraightVelocityUpdater) correctVelocitiesSampling(_ Spheres, _ identifier, _ *systemParameters) {
	return
}

type ReflectiveVelocityUpdater struct{}

func (u *ReflectiveVelocityUpdater) initializeVelocities(spheres Spheres, sys *systemParameters) identifier {
	actId := identifier(rand.Intn(len(spheres)))
	randomPhi := -math.Pi + 2.0*math.Pi*rand.Float64()
	velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomPhi))
	velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomPhi))
	spheres[actId].velocity = []*big.Float{velocityX, velocityY}
	return actId
}

func (u *ReflectiveVelocityUpdater) updateVelocityEvent(activeSphere, targetSphere *Sphere, sys *systemParameters) {
	separation := separationVector(activeSphere.position, targetSphere.position, sys)
	separationNorm := new(big.Float)
	separationNorm.Sqrt(normSquared(separation))
	for index := range separation {
		separation[index].Quo(separation[index], separationNorm)
	}
	velocityDotSeparation := dot(separation, activeSphere.velocity)
	velocityX := new(big.Float)
	velocityY := new(big.Float)
	two := new(big.Float).SetPrec(sys.prec).SetInt64(2)
	velocityX.Mul(separation[0], velocityDotSeparation)
	velocityX.Mul(velocityX, two)
	velocityX.Sub(velocityX, activeSphere.velocity[0])
	velocityY.Mul(separation[1], velocityDotSeparation)
	velocityY.Mul(velocityY, two)
	velocityY.Sub(velocityY, activeSphere.velocity[1])
	targetSphere.velocity = []*big.Float{velocityX, velocityY}
	velocityNorm := normSquared(targetSphere.velocity)
	velocityNorm.Sqrt(velocityNorm)
	for index := range targetSphere.velocity {
		targetSphere.velocity[index].Quo(targetSphere.velocity[index], velocityNorm)
	}
	activeSphere.velocity = nil
}

func (u *ReflectiveVelocityUpdater) updateVelocitiesEndOfChain(spheres Spheres, actId identifier,
	sys *systemParameters) identifier {

	newActId := identifier(rand.Intn(len(spheres)))
	randomPhi := -math.Pi + 2.0*math.Pi*rand.Float64()
	velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomPhi))
	velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomPhi))
	spheres[newActId].velocity = []*big.Float{velocityX, velocityY}
	if newActId != actId {
		spheres[actId].velocity = nil
	}
	return newActId
}

func (u *ReflectiveVelocityUpdater) correctVelocitiesSampling(_ Spheres, _ identifier, _ *systemParameters) {
	return
}

type ForwardVelocityUpdater struct{}

func (u *ForwardVelocityUpdater) initializeVelocities(spheres Spheres, sys *systemParameters) identifier {
	actId := identifier(rand.Intn(len(spheres)))
	randomPhi := -math.Pi + 2.0*math.Pi*rand.Float64()
	velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomPhi))
	velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomPhi))
	spheres[actId].velocity = []*big.Float{velocityX, velocityY}
	return actId
}

func (u *ForwardVelocityUpdater) updateVelocityEvent(activeSphere, targetSphere *Sphere, sys *systemParameters) {
	// After normalization, n is the unit vector in the direction of the gradient of the sphere potential.
	n := separationVector(activeSphere.position, targetSphere.position, sys)
	nNorm := new(big.Float)
	nNorm.Sqrt(normSquared(n))
	for index := range n {
		n[index].Quo(n[index], nNorm)
	}
	nPerpendicularX := new(big.Float)
	nPerpendicularY := new(big.Float)
	nPerpendicularX.Neg(n[1])
	nPerpendicularY.Copy(n[0])
	// Unit vector perpendicular to n.
	nPerpendicular := []*big.Float{nPerpendicularX, nPerpendicularY}
	// We randomly sample the parts of the new active particle's velocity that are parallel and perpendicular to n.
	orthogonalPart := rand.Float64()
	parallelPart := math.Sqrt(1.0 - orthogonalPart*orthogonalPart)
	// The parallel part of the new active particle's velocity should have the same sign as before which corresponds to
	// inverting the sign of the parallel part if the old active particle would have stayed active.
	if dot(activeSphere.velocity, n).Sign() < 0 {
		parallelPart = -parallelPart
	}
	// The sign of the orthogonal part of the new active particle's velocity should be inverted which corresponds to
	// keeping the sign of the orthogonal part if the old active particle would have stayed active.
	if dot(activeSphere.velocity, nPerpendicular).Sign() >= 0 {
		orthogonalPart = -orthogonalPart
	}
	velocityX := new(big.Float)
	velocityY := new(big.Float)
	parallel := new(big.Float).SetPrec(sys.prec).SetFloat64(parallelPart)
	orthogonal := new(big.Float).SetPrec(sys.prec).SetFloat64(orthogonalPart)
	t := new(big.Float)
	velocityX.Mul(parallel, n[0])
	t.Mul(orthogonal, nPerpendicularX)
	velocityX.Add(velocityX, t)
	velocityY.Mul(parallel, n[1])
	t.Mul(orthogonal, nPerpendicularY)
	velocityY.Add(velocityY, t)
	targetSphere.velocity = []*big.Float{velocityX, velocityY}
	velocityNorm := normSquared(targetSphere.velocity)
	velocityNorm.Sqrt(velocityNorm)
	for index := range targetSphere.velocity {
		targetSphere.velocity[index].Quo(targetSphere.velocity[index], velocityNorm)
	}
	activeSphere.velocity = nil
}

func (u *ForwardVelocityUpdater) updateVelocitiesEndOfChain(spheres Spheres, actId identifier,
	sys *systemParameters) identifier {

	newActId := identifier(rand.Intn(len(spheres)))
	randomPhi := -math.Pi + 2.0*math.Pi*rand.Float64()
	velocityX := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Cos(randomPhi))
	velocityY := new(big.Float).SetPrec(sys.prec).SetFloat64(math.Sin(randomPhi))
	spheres[newActId].velocity = []*big.Float{velocityX, velocityY}
	if newActId != actId {
		spheres[actId].velocity = nil
	}
	return newActId
}

func (u *ForwardVelocityUpdater) correctVelocitiesSampling(_ Spheres, _ identifier, _ *systemParameters) {
	return
}

type NewtonianVelocityUpdater struct {
	StandardDeviation float64
}

func (u *NewtonianVelocityUpdater) sampleVelocitiesMaxwellBoltzmann(spheres Spheres, sys *systemParameters) {
	var sphereIndex, v int
	var sphere *Sphere
	var randEntry float64
	for sphereIndex = range spheres {
		sphere = &spheres[sphereIndex]
		sphere.velocity = make([]*big.Float, 2)
		for v = range sphere.velocity {
			randEntry = rand.NormFloat64() * u.StandardDeviation
			sphere.velocity[v] = new(big.Float).SetPrec(sys.prec).SetFloat64(randEntry)
		}
	}
	u.correctVelocitiesMaxwellBoltzmann(spheres, sys)
}

func (u *NewtonianVelocityUpdater) correctVelocitiesMaxwellBoltzmann(spheres Spheres, sys *systemParameters) {
	var sphereIndex, v int
	var sphere *Sphere
	meanVelocity := []*big.Float{new(big.Float).SetPrec(sys.prec), new(big.Float).SetPrec(sys.prec)}
	t := new(big.Float)
	nSpheresFloat := new(big.Float).SetPrec(sys.prec).SetInt64(int64(len(spheres)))
	for sphereIndex = range spheres {
		sphere = &spheres[sphereIndex]
		for v = range sphere.velocity {
			meanVelocity[v].Add(meanVelocity[v], t.Quo(sphere.velocity[v], nSpheresFloat))
		}
	}

	meanVelocitySquared := new(big.Float).SetPrec(sys.prec)
	// Correct overall momentum to zero.
	for sphereIndex = range spheres {
		sphere = &spheres[sphereIndex]
		for v = range sphere.velocity {
			sphere.velocity[v].Sub(sphere.velocity[v], meanVelocity[v])
			t.Mul(sphere.velocity[v], sphere.velocity[v])
			meanVelocitySquared.Add(meanVelocitySquared, t.Quo(t, nSpheresFloat))
		}
	}

	// Correct root mean square velocity to 1.
	rootMeanSquareVelocity := new(big.Float)
	rootMeanSquareVelocity.Sqrt(meanVelocitySquared)
	for sphereIndex = range spheres {
		sphere = &spheres[sphereIndex]
		for v = range sphere.velocity {
			sphere.velocity[v].Quo(sphere.velocity[v], rootMeanSquareVelocity)
		}
	}
}

func (u *NewtonianVelocityUpdater) initializeVelocities(spheres Spheres, sys *systemParameters) identifier {
	u.sampleVelocitiesMaxwellBoltzmann(spheres, sys)
	return identifier(rand.Intn(len(spheres)))
}

func (u *NewtonianVelocityUpdater) updateVelocityEvent(activeSphere, targetSphere *Sphere, sys *systemParameters) {
	separation := separationVector(activeSphere.position, targetSphere.position, sys)
	separationSquared := normSquared(separation)
	velocityDifferenceDotSeparation := new(big.Float).SetPrec(sys.prec)
	t := new(big.Float)
	var vIndex int
	for vIndex = range activeSphere.velocity {
		t.Sub(targetSphere.velocity[vIndex], activeSphere.velocity[vIndex])
		velocityDifferenceDotSeparation.Add(velocityDifferenceDotSeparation, t.Mul(t, separation[vIndex]))
	}
	t2 := new(big.Float)
	t2.Quo(velocityDifferenceDotSeparation, separationSquared)
	for vIndex = range activeSphere.velocity {
		t.Mul(separation[vIndex], t2)
		activeSphere.velocity[vIndex].Add(activeSphere.velocity[vIndex], t)
		targetSphere.velocity[vIndex].Sub(targetSphere.velocity[vIndex], t)
	}
}

func (u *NewtonianVelocityUpdater) updateVelocitiesEndOfChain(spheres Spheres, _ identifier,
	sys *systemParameters) identifier {

	u.sampleVelocitiesMaxwellBoltzmann(spheres, sys)
	return identifier(rand.Intn(len(spheres)))
}

func (u *NewtonianVelocityUpdater) correctVelocitiesSampling(spheres Spheres, _ identifier, sys *systemParameters) {
	u.correctVelocitiesMaxwellBoltzmann(spheres, sys)
}

func runChainCells(spheres Spheres, actId identifier, chainTime *big.Float, sys *systemParameters, c *Cells,
	velocityUpdater VelocityUpdater, numberEvents *uint64, samplingTime *big.Float,
	remainingSamplingTime *big.Float, file File, maximumNumberEvents uint64) identifier {

	activeSphere := &spheres[actId]
	activeCell := c.positionToCell(activeSphere.position)
	timeToGo := new(big.Float)
	timeToGo.Copy(chainTime)
	soonestEvent := event{new(big.Float).SetPrec(sys.prec).SetInt64(-1), identifier(-1)}
	var nextCell cellId
	var time, entry *big.Float
	var positionIndex int
	var targetId identifier
	t := new(big.Float)
	for {
		// Compute cell boundary event.
		time, nextCell = c.timeToCellBoundary(activeSphere, activeCell, sys)
		soonestEvent.time.Copy(time)
		soonestEvent.nextId = actId

		// Compute sphere events.
		for _, targetId = range c.occupantsNeighborCells(activeCell) {
			if targetId != actId {
				time = timeToSphereEvent(activeSphere, &spheres[targetId], sys)
				if time.Sign() < 0 {
					time = new(big.Float).SetPrec(sys.prec)
				}
				if time.Cmp(soonestEvent.time) < 0 {
					soonestEvent.time.Copy(time)
					soonestEvent.nextId = targetId
				}
			}
		}

		// Sampling event.
		if remainingSamplingTime.Cmp(timeToGo) <= 0 && remainingSamplingTime.Cmp(soonestEvent.time) <= 0 {
			for positionIndex, entry = range activeSphere.position {
				entry.Add(entry, t.Mul(remainingSamplingTime, activeSphere.velocity[positionIndex]))
				if entry.Cmp(sys.systemLength) >= 0 {
					entry.Sub(entry, sys.systemLength)
				} else if entry.Sign() < 0 {
					entry.Add(entry, sys.systemLength)
				}
			}
			file.Write(spheres, sys, *numberEvents)
			timeToGo.Sub(timeToGo, remainingSamplingTime)
			remainingSamplingTime.Copy(samplingTime)
			velocityUpdater.correctVelocitiesSampling(spheres, actId, sys)
			if file.Finished() {
				return actId
			}
			continue
		}

		// End of chain.
		if timeToGo.Cmp(soonestEvent.time) <= 0 {
			*numberEvents += 1
			for positionIndex, entry = range activeSphere.position {
				entry.Add(entry, t.Mul(timeToGo, activeSphere.velocity[positionIndex]))
				if entry.Cmp(sys.systemLength) >= 0 {
					entry.Sub(entry, sys.systemLength)
				} else if entry.Sign() < 0 {
					entry.Add(entry, sys.systemLength)
				}
			}
			remainingSamplingTime.Sub(remainingSamplingTime, timeToGo)
			if *numberEvents%10000 == 0 {
				checkState(spheres, sys, false)
			}
			if *numberEvents > maximumNumberEvents {
				return -1
			}
			break
		}

		// Real or cell-boundary event.
		for positionIndex, entry = range activeSphere.position {
			entry.Add(entry, t.Mul(soonestEvent.time, activeSphere.velocity[positionIndex]))
			if entry.Cmp(sys.systemLength) >= 0 {
				entry.Sub(entry, sys.systemLength)
			} else if entry.Sign() < 0 {
				entry.Add(entry, sys.systemLength)
			}
		}
		timeToGo.Sub(timeToGo, soonestEvent.time)
		remainingSamplingTime.Sub(remainingSamplingTime, soonestEvent.time)
		if actId == soonestEvent.nextId {
			// Cell-boundary event. Update cell-occupancy system and change cell of active sphere.
			c.updateCell(actId, activeCell, nextCell)
			activeCell = nextCell
		} else {
			// Real event, active sphere changed.
			*numberEvents += 1
			velocityUpdater.updateVelocityEvent(activeSphere, &spheres[soonestEvent.nextId], sys)
			actId = soonestEvent.nextId
			activeSphere = &spheres[actId]
			activeCell = c.positionToCell(activeSphere.position)
			if *numberEvents%10000 == 0 {
				checkState(spheres, sys, false)
			}
			if *numberEvents > maximumNumberEvents {
				return -1
			}
		}
	}
	return actId
}
