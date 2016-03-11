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
#ifndef DG_H_
#define DG_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "SmbData.h"
#include "math/Field.h"
#include "solver/SpatialDiscretization.h"
#include "solver/dgtd/core/Comm.h"
#include "solver/dgtd/core/BCGroup.h"
#include "options/OptionsSolverDGTD.h"
#include "../dg/sources/DGPlaneWave.h"
#include "../dg/sources/DGDipole.h"
#include "../dg/sources/DGWaveguidePortRectangular.h"
#include "../dg/dispersive/DGSIBC.h"
#include "../dg/dispersive/DGDispersiveVolumic.h"
//#include "dispersives/DGPMLUniaxial.h"
//#include "dispersives/DGPMLMultiaxial.h"

#define SOLVER_DEDUPLICATE_OPERATORS

struct lexCompareMat {
    static const UInt np = ((ORDER_N+1) * (ORDER_N+2) * (ORDER_N+3) / 6);
    bool
    operator() (
            const StaMatrix<Real,np,np>& lhs,
            const StaMatrix<Real,np,np>& rhs) const {
        static const Real tolerance = 1e-12;
        for (UInt i = 0; i < (np*np); i++) {
            if (std::abs(lhs.val(i) - rhs.val(i)) > tolerance) {
                if (lhs.val(i) < rhs.val(i)) {
                    return true;
                }
                if (lhs.val(i) > rhs.val(i)) {
                    return false;
                }
            }
        }
        return false;
    }
};

class DG : public Ordering, public SpatialDiscretization {
    friend class Exporter;
    friend class IntegratorLSERK;
    friend class IntegratorLF2;
    friend class IntegratorLF2Full;
    friend class IntegratorVerlet;
public:
    const static UInt N = ORDER_N;
    const static UInt faces = 4;
    const static UInt np = (N+1) * (N+2) * (N+3) / 6;
    const static UInt np2 = np * 2;
    const static UInt nfp = (N+1) * (N+2) / 2;
    const static UInt npnfp = np * nfp;
    const static UInt npnp = np * np;
    const static UInt nfpfaces = nfp * faces;
    DG();
    virtual ~DG();
    virtual void setFieldsToRandom();
    virtual void setFieldsToGaussian(
            const CellGroup& cells,
            const Real amplitude,
            CVecR3& polarization,
            const CVecR3& gaussCenter,
            const Real gaussWidth);
    virtual void setFieldsToHarmonics(
            const CellGroup& cells,
            const CVecI3& harmonics,
            CVecR3& polarization);
    void setFieldsAndTimeFromResumeFile();
    virtual UInt getFieldDOFs() = 0;
    const FieldR3* getElectric() const;
    const FieldR3* getMagnetic() const;
    virtual UInt getGlobalFieldPosOfVertex(
            pair<const ElemR*, UInt> vertex) const;
    virtual vector<UInt> getGlobalFieldPosOfFace(
            Face boundary) const;
    virtual vector<UInt> getGlobalFieldPosOfVolume(
            const ElemId volId) const;
    virtual void printInfo() const = 0;

protected:
    FieldR3 E, H;
    FieldR3 resE, resH;
    FieldR3 nAdm, nImp, rnAdm, rnImp, cnAdm, cnImp;
    Real *oneOverEps, *oneOverMu;
    Comm* comm;
    UInt nK;
    Real upwinding;
    Real LIFT[faces*npnfp]; //!< Flux gatherer operator. dim = matrix(np x (4*nfp))

#ifdef SOLVER_DEDUPLICATE_OPERATORS
    set<StaMatrix<Real,np,np>, lexCompareMat> CList;
#else
    StaMatrix<Real,np,np>* CList;
#endif
    void init(
            const OptionsSolverDGTD& options,
            const PMGroup& pm,
            const CellGroup& cells,
            Comm* comm_);
    virtual void addFluxesToRHS(
            const UInt e1, const UInt e2,
            const Real localTime,
            const Real minDT);
    virtual void addFluxesToRHSElectric(const UInt e1, const UInt e2) = 0;
    virtual void addFluxesToRHSMagnetic(const UInt e1, const UInt e2) = 0;
    virtual void addFluxesToRHSElectric(
            const UInt e1, const UInt e2, const bool useResForUpw) = 0;
    virtual void addFluxesToRHSMagnetic(
            const UInt e1, const UInt e2, const bool useResForUpw) = 0;
    virtual void addRHSToFieldsElectric(
            const UInt e1,
            const UInt e2,
            const Real rkdt) = 0;
    virtual void addRHSToFieldsMagnetic(
            const UInt e1,
            const UInt e2,
            const Real rkdt) = 0;
    virtual void addRHSToResidueElectric(
            const UInt e1, const UInt e2, const Real rkdt) = 0;
    virtual void addRHSToResidueMagnetic(
            const UInt e1, const UInt e2, const Real rkdt) = 0;
    void buildCMatrices(const CellGroup& cells);
    virtual void buildMaterials(
            const CellGroup& cells,
            const OptionsSolverDGTD& arg) = 0;
    virtual void  computeRHS(
            const UInt e1,
            const UInt e2,
            const Real localtime,
            const Real rkdt) = 0;
    virtual void computeRHSElectric(
            const UInt e1,
            const UInt e2,
            const Real localtime,
            const Real minDT) = 0;
    virtual void computeRHSMagnetic(
            const UInt e1,
            const UInt e2,
            const Real localtime,
            const Real minDT) = 0;
    virtual UInt getIndexOfElement(const UInt e) const = 0;
    virtual const FieldR3& getRHSElectric() const = 0;
    virtual const FieldR3& getRHSMagnetic() const = 0;
    virtual void updateFieldsWithRes(
            const UInt e1,
            const UInt e2,
            const Real rkb) = 0;
    virtual void updateFieldsWithResBase(
            const UInt e1,
            const UInt e2,
            const Real rkb);
    virtual void computeCurlsInRHS(const UInt e1, const UInt e2);
    virtual void computeCurlsInRHSElectric(const UInt e1, const UInt e2) = 0;
    virtual void computeCurlsInRHSMagnetic(const UInt e1, const UInt e2) = 0;
    virtual void computeJumps(
            const UInt e1,
            const UInt e2,
            const Real localTime,
            const Real minDT) = 0;
    void copyFieldsInResidues(
            const UInt e1,
            const UInt e2);
    virtual void copyJumpsToResidueJumps(
            const UInt e1,
            const UInt e2) = 0;
    virtual void computePolarizationCurrentsRHS(
            const UInt e1, const UInt e2) = 0;
    virtual void computePolarizationCurrentsRHSElectric(
            const UInt e1, const UInt e2) = 0;
    virtual void computePolarizationCurrentsRHSMagnetic(
            const UInt e1, const UInt e2) = 0;
    virtual void LTSSaveFieldsAndResidues(
            const UInt fKSave,
            const UInt lKSave) = 0;
    virtual void LTSLoadFieldsAndResidues(
            const UInt fKSave,
            const UInt lKSave) = 0;
    virtual void addRHSToRes(
            const UInt e1,
            const UInt e2,
            const Real rka,
            const Real dt) = 0;
    void swapResiduesAndFields(
            const UInt e1,
            const UInt e2);
    void buildFluxScalingFactors(const CellGroup& cells, const Connectivities& map);
    void buildFieldScalingFactors(const CellGroup& cells);
    virtual void buildScalingFactors(
            const CellGroup& cells,
            const Connectivities& map);
    void buildLIFT();
    virtual void assignMatrices(const CellGroup& cells) = 0;
    void allocateFieldsAndRes();
};
#endif /* SOLVER_H_ */

