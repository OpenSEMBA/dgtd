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

#include "source/port/WaveguideRectangular.h"
#include "Waveport.h"

namespace SEMBA {
namespace Cudg3d {

//class DGWaveportRectangular : public DGWaveport, public Source::Port::WaveguideRectangular {
//public:
//    DGWaveportRectangular(
//            const PortWaveguide& pw,
//            const Connectivities& map,
//            FieldR3& dE, FieldR3& dH,
//            const Math::Int vmapM[faces][nfp]);
//    virtual ~DGWaveportRectangular();
//    void computeExcitation(
//            const Math::Real intTime,
//            const Math::Real minDT);
//    void printInfo() const;
//private:
//    Math::Real width, height;
//    Source::Port::Waveguide::ExcitationMode excitationMode;
//    Math::Real kcm;
//    Math::Real intrinsicImpedance;
//    Math::Real gammaMSum;
//    void computeExcitationField(
//            FieldR3& EInc,
//            FieldR3& HInc,
//            const Math::CVecR3* pos,
//            const size_t nE,
//            const Math::Real intTime,
//            const Math::Real minDT);
//};

}

#endif /* SOLVERWAVEPORTRECTANGULAR_H_ */
