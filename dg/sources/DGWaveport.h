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
 * SolverWaveport.h
 *
 *  Created on: Aug 26, 2013
 *      Author: luis
 */

#ifndef SOLVERWAVEPORT_H_
#define SOLVERWAVEPORT_H_

#include "../../dg/sources/DGSource.h"

class DGWaveport : public DGSource {
public:
    DGWaveport();
    virtual	~DGWaveport();
protected:
    //	Real
    //     getNumericalGammaMGauss(
    //      const Real time,
    //      const Real minDT,
    //      const Real amplitude,
    //      const Real delay,
    //      const Real spread,
    //      const Real kcm) const;
    bool checkNormalsAreEqual(
            const vector<pair<UInt,UInt> >& elemFace) const;
protected:
    CVecR3* posTF;
    CVecR3* posTFNB;
    CVecR3* posSF;
private:
    Real *gauss, *hm;
    Real getHm(
            const Real time,
            const Real kcm) const;
};

#endif /* SOLVERWAVEPORT_H_ */
