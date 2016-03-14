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
 * SolverLSERK.h
 *
 *  Created on: Nov 30, 2012
 *      Author: luis
 */
#ifndef INTEGRATORLSERK_H_
#define INTEGRATORLSERK_H_

#include "../../dgtd/integrator/Integrator.h"

class IntegratorLSERK : public Integrator {
    friend class DG;
public:
    IntegratorLSERK();
    virtual ~IntegratorLSERK();
    IntegratorLSERK(
            const Mesh::Volume& mesh,
            const PMGroup& pmGroup,
            const OptionsSolverDGTD* arg);
    void timeIntegrate(const Real time) const;
protected:
    size_t getNumOfIterationsPerBigTimeStep(
            const size_t e) const;
private:
    bool useMaxStageSizeForLTS;
    static const size_t nStages = 5;
    static const Real rka[nStages];
    static const Real rkb[nStages];
    static const Real rkc[nStages];
    Real stageSize[nStages];
    Real getMaxStageSize() const;
    Real getMaxTimeRatio() const;
    void buildRKConstants();
    void LTSTimeIntegration(
            const Real time,
            const Real localTime,
            const Real localdt,
            const size_t tier) const;
    void updateResiduals(
            const size_t e1,
            const size_t e2,
            const Real rka,
            const Real localTime,
            const Real localdt) const;
    size_t getNStages() const;
    Real getRKA(const size_t s) const;
    Real getRKB(const size_t s) const;
    Real getRKC(const size_t s) const;
    Real getStageSize(const size_t s) const;
};
#endif /* SOLVERLSERK_H_ */
