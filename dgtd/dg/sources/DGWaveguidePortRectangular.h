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
 * SolverWaveportRectangular.h
 *
 *  Created on: Aug 26, 2013
 *      Author: luis
 */

#ifndef SOLVERWAVEPORTRECTANGULAR_H_
#define SOLVERWAVEPORTRECTANGULAR_H_

#include "../../dg/sources/DGWaveport.h"
#include "sources/ports/PortWaveguideRectangular.h"

class DGWaveportRectangular : public DGWaveport, public PortWaveguideRectangular {
public:
    DGWaveportRectangular(
            const PortWaveguide& pw,
            const Connectivities& map,
            FieldR3& dE, FieldR3& dH,
            const Int vmapM[faces][nfp]);
    virtual ~DGWaveportRectangular();
    void computeExcitation(
            const Real intTime,
            const Real minDT);
    void printInfo() const;
private:
    Real width, height;
    PortWaveguide::ExcitationMode excitationMode;
    Real kcm;
    Real intrinsicImpedance;
    Real gammaMSum;
    void computeExcitationField(
            FieldR3& EInc,
            FieldR3& HInc,
            const CVecR3* pos,
            const UInt nE,
            const Real intTime,
            const Real minDT);
};

#endif /* SOLVERWAVEPORTRECTANGULAR_H_ */
