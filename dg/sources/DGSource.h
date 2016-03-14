// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
/*
 * SolverSource.h
 *
 *  Created on: Sep 2, 2013
 *      Author: luis
 */

#ifndef DG_SOLVER_SOURCE_H_
#define DG_SOLVER_SOURCE_H_

#include <utility>
#include <vector>

using namespace std;

#include "sources/EMSource.h"
#include "communications/Comm.h"
#include "boundaryConditions/Group.h"

class DGSource {
public:
    const static size_t N = ORDER_N;
    const static size_t nfp = (N+1) * (N+2) / 2;
    const static size_t faces = 4;
    typedef enum {
        totalField,
        scatteredField,
        totalFieldNotBacked
    } BackingType;
    DGSource();
    virtual ~DGSource();
    void addJumps(
            const size_t e1,
            const size_t e2);
    virtual void computeExcitation(
            const Real intTime,
            const Real minDT) = 0;
    virtual void printInfo() const = 0;
protected:
    const static size_t np = (N+1) * (N+2) * (N+3) / 6;
    const static size_t np2 = np * 2;
    const static size_t npnfp = np * nfp;
    const static size_t npnp = np * np;
    const static size_t nfpfaces = nfp * faces;
    // Excitation fields.
    FieldR3 ETInc, ESInc, EIncNB;
    FieldR3 HTInc, HSInc, HIncNB;

    vector<size_t> ETFe, ESFe, ETFNBe;
    // Excitation total field jumps pointers.
    Real **dExT, **dEyT, **dEzT;
    Real **dHxT, **dHyT, **dHzT;
    // Excitation scattered field jumps pointers.
    Real **dExS, **dEyS, **dEzS;
    Real **dHxS, **dHyS, **dHzS;
    // Excitation total field not backed jumps.
    Real **dExTNB, **dEyTNB, **dEzTNB;
    Real **dHxTNB, **dHyTNB, **dHzTNB;
    void initSource(
            const BCGroup& bc,
            const Connectivities& map,
            const CellGroup& cells,
            FieldR3& dE,
            FieldR3& dH,
            const Int vmapM[faces][nfp]);
    CVecR3* initPositions(
            const vector<pair<size_t, size_t> >& elemFace,
            const CellGroup& cells) const;
    vector<pair<size_t, size_t>> getTotalFieldElemFaces(
            const BCGroup& bc,
            const Connectivities& map,
            const CellGroup& cells) const;
    vector<pair<size_t, size_t>> getScattFieldElemFaces(
            const BCGroup& bc,
            const Connectivities& map,
            const CellGroup& cells) const;
    vector<pair<size_t, size_t>> getTotalNotBackedFieldElemFaces(
            const BCGroup& bc,
            const Connectivities& map,
            const CellGroup& cells) const;
};

#endif /* SOLVERSOURCE_H_ */
