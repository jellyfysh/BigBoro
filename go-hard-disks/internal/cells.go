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

type Cells struct {
	occupants      [][]identifier
	cellsPerSide   int64
	cellSideLength *big.Float
}

type cellId [2]int64

func (cells *Cells) cellIndex(cellId cellId) int64 {
	return cellId[0] + cellId[1]*cells.cellsPerSide
}

func (cells *Cells) positionToCell(position []*big.Float) cellId {
	q := new(big.Float)
	idx, _ := q.Quo(position[0], cells.cellSideLength).Int64()
	idy, _ := q.Quo(position[1], cells.cellSideLength).Int64()
	return cellId{idx, idy}
}

func CreateCells(spheres Spheres, cellsPerSide int, sys *systemParameters) *Cells {
	if !(cellsPerSide > 2) {
		panic("Cell-occupancy system only allows for at least 3 cells per side.")
	}
	t := new(big.Float).SetPrec(sys.prec).SetInt64(int64(cellsPerSide))
	if t.Quo(sys.systemLength, t).Cmp(sys.twoSigma) < 0 {
		panic("Too many cells per side.")
	}
	cellSideLength := new(big.Float)
	cellSideLength.Quo(sys.systemLength, new(big.Float).SetPrec(sys.prec).SetInt64(int64(cellsPerSide)))
	emptyOccupants := make([][]identifier, cellsPerSide*cellsPerSide)
	for i := range emptyOccupants {
		emptyOccupants[i] = make([]identifier, 0, 8)
	}
	cells := Cells{
		occupants:      emptyOccupants,
		cellsPerSide:   int64(cellsPerSide),
		cellSideLength: cellSideLength,
	}
	var cellIndex int64
	for sphereIndex := range spheres {
		cellIndex = cells.cellIndex(cells.positionToCell(spheres[sphereIndex].position))
		cells.occupants[cellIndex] = append(
			cells.occupants[cellIndex], identifier(sphereIndex))
	}
	return &cells
}

func (cells *Cells) updateCell(activeIdentifier identifier, oldCell, newCell cellId) {
	removed := false
	oldCellIndex := cells.cellIndex(oldCell)
	newCellIndex := cells.cellIndex(newCell)
	for index, id := range cells.occupants[oldCellIndex] {
		if id == activeIdentifier {
			removed = true
			cells.occupants[oldCellIndex] = append(
				cells.occupants[oldCellIndex][:index], cells.occupants[oldCellIndex][index+1:]...)
			break
		}
	}

	if !removed {
		var realCellId cellId
		found := false
		var cellIndex int64
		for cellIndex = 0; cellIndex < cells.cellsPerSide*cells.cellsPerSide; cellIndex++ {
			for _, id := range cells.occupants[cellIndex] {
				if id == activeIdentifier {
					realCellId = cellId{cellIndex % cells.cellsPerSide, cellIndex / cells.cellsPerSide}
					found = true
					break
				}
			}
			if found {
				break
			}
		}
		panic(fmt.Sprintf("Active identifier not found in given old cell %v but in cell %v (new cell %v).",
			oldCell, realCellId, newCell))
	}

	cells.occupants[newCellIndex] = append(cells.occupants[newCellIndex], activeIdentifier)
}

func (cells *Cells) occupantsNeighborCells(activeCell cellId) []identifier {
	occ := make([]identifier, 0, 9)
	var id identifier
	var nIndex, j int64
	var nxIndices, nyIndices [3]int64

	switch activeCell[0] {
	case 0:
		nxIndices = [3]int64{cells.cellsPerSide - 1, 0, 1}
	case cells.cellsPerSide - 1:
		nxIndices = [3]int64{cells.cellsPerSide - 2, cells.cellsPerSide - 1, 0}
	default:
		nxIndices = [3]int64{activeCell[0] - 1, activeCell[0], activeCell[0] + 1}
	}

	switch activeCell[1] {
	case 0:
		nyIndices = [3]int64{cells.cellsPerSide - 1, 0, 1}
	case cells.cellsPerSide - 1:
		nyIndices = [3]int64{cells.cellsPerSide - 2, cells.cellsPerSide - 1, 0}
	default:
		nyIndices = [3]int64{activeCell[1] - 1, activeCell[1], activeCell[1] + 1}
	}

	for i := 0; i < 3; i++ {
		for j = 0; j < 3; j++ {
			nIndex = nxIndices[i] + nyIndices[j]*cells.cellsPerSide
			for _, id = range cells.occupants[nIndex] {
				occ = append(occ, id)
			}
		}
	}
	return occ
}

func (cells *Cells) timeToCellBoundary(activeSphere *Sphere, activeCell cellId,
	sys *systemParameters) (*big.Float, cellId) {

	boundary := new(big.Float)
	separation := new(big.Float)
	time := new(big.Float)
	var shortestDir int
	// Positive infinity.
	shortestTime := new(big.Float).SetPrec(sys.prec).SetInf(false)
	nextCell := activeCell
	t := new(big.Float).SetPrec(sys.prec)

	for i, v := range activeSphere.velocity {
		if v.Sign() != 0 {
			if v.Sign() > 0 {
				boundary.Mul(t.SetInt64(activeCell[i]+1), cells.cellSideLength)
				// Position 0.0 might be stored as L due to periodic boundary conditions.
				if activeSphere.position[i].Cmp(boundary) > 0 {
					boundary.Add(boundary, sys.systemLength)
				}
			} else {
				boundary.Mul(t.SetInt64(activeCell[i]), cells.cellSideLength)
				// Position L might be stored as 0.0 due to periodic boundary conditions.
				if activeSphere.position[i].Cmp(boundary) < 0 {
					boundary.Sub(boundary, sys.systemLength)
				}
			}
			separation.Sub(boundary, activeSphere.position[i])
			time.Quo(separation, v)
			if time.Cmp(shortestTime) < 0 {
				shortestTime.Copy(time)
				shortestDir = i
			}
		}
	}

	activeCellDir := activeCell[shortestDir]
	if activeSphere.velocity[shortestDir].Sign() > 0 {
		if activeCellDir < cells.cellsPerSide-1 {
			nextCell[shortestDir] += 1
		} else {
			nextCell[shortestDir] = 0
		}
	} else {
		if activeCellDir > 0 {
			nextCell[shortestDir] -= 1
		} else {
			nextCell[shortestDir] = cells.cellsPerSide - 1
		}
	}
	return shortestTime, nextCell
}
