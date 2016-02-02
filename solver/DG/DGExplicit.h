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
 * SolverExplicit.h
 *
 *  Created on: Nov 30, 2012
 *      Author: luis
 */

#ifndef SOLVEREXPLICIT_H_

#include "DG.h"
#include "DGCurvedFace.h"
#include "SolverMath.h"

class DGExplicit : public DG {
    friend class IntegratorLSERK;
    friend class IntegratorLF2;
    friend class IntegratorLF2Full;
    friend class IntegratorVerlet;
public:
    DGExplicit(
            const MeshVolume& mesh,
            const PMGroup& pMGroup,
            const EMSourceGroup& emsources,
            const OptionsSolverDGTD& options,
            Comm* comm);
    virtual ~DGExplicit();
    UInt getFieldDOFs();
    const FieldR3& getRHSElectric() const;
    const FieldR3& getRHSMagnetic() const;
    void printInfo() const;
protected:
    void computeRHS(
            const UInt e1,
            const UInt e2,
            const Real localtime,
            const Real rkdt);
    void computeRHSElectric(
            const UInt e1,
            const UInt e2,
            const Real localtime,
            const Real minDT);
    void computeRHSMagnetic(
            const UInt e1,
            const UInt e2,
            const Real localtime,
            const Real minDT);
    void computeCurlsInRHSElectric(const UInt e1, const UInt e2);
    void computeCurlsInRHSMagnetic(const UInt e1, const UInt e2);
    void computeJumps(
            const UInt e1,
            const UInt e2,
            const Real localTime,
            const Real minDT);
    void copyJumpsToResidueJumps(
            const UInt e1,
            const UInt e2);
    void addFluxesToRHSElectric(
            const UInt e1, const UInt e2);
    void addFluxesToRHSMagnetic(
            const UInt e1, const UInt e2);
    void addFluxesToRHSElectric(
            const UInt e1, const UInt e2, const bool useResForUpw);
    void addFluxesToRHSMagnetic(
            const UInt e1, const UInt e2, const bool useResForUpw);
    void addStraightFluxesToRHSElectric(const UInt e1, const UInt e2,
            const bool useResForUpw);
    void addStraightFluxesToRHSMagnetic(const UInt e1, const UInt e2,
            const bool useResForUpw);
    void addCurvedFluxesToRHSElectric(const UInt e1, const UInt e2,
            const bool useResForUpw);
    void addCurvedFluxesToRHSMagnetic(const UInt e1, const UInt e2,
            const bool useResForUpw);
    void computePolarizationCurrentsRHS(
            const UInt e1, const UInt e2);
    void computePolarizationCurrentsRHSElectric(
            const UInt e1, const UInt e2);
    void computePolarizationCurrentsRHSMagnetic(
            const UInt e1, const UInt e2);
    void addRHSToFieldsElectric(
            const UInt e1,
            const UInt e2,
            const Real rkdt);
    void addRHSToFieldsMagnetic(
            const UInt e1,
            const UInt e2,
            const Real rkdt);
    UInt getIndexOfElement(const UInt e) const;
    void addRHSToResidueElectric(const UInt e1, const UInt e2,
            const Real rkdt);
    void addRHSToResidueMagnetic(const UInt e1, const UInt e2,
            const Real rkdt);
    void addRHSToRes(
            const UInt e1,
            const UInt e2,
            const Real rka,
            const Real dt);
    void updateFieldsWithRes(
            const UInt e1,
            const UInt e2,
            const Real rkb);
    void LTSSaveFieldsAndResidues(
            const UInt fKSave,
            const UInt lKSave);
    void LTSLoadFieldsAndResidues(
            const UInt fKSave,
            const UInt lKSave);
private:
    // - Maps.
    Int vmapM[faces][nfp];
    Int vmapP[16][nfp];
    Int ***map_;
    // Pointers to neighbour fields. dim = (nK, 4).
    Real ***ExP, ***EyP, ***EzP, ***HxP, ***HyP, ***HzP;
    // Curved faces stuff ---------------------------------------------
    UInt nCurvedFaces;
    DGCurvedFace *curveFace;
    const Real **Cx, **Cy, **Cz; // Pointers to C. dim = (nK)
    // Fields and residuals: dim = (np,nK)
    FieldR3 rhsE, rhsH;
    FieldR3 savedResE, savedResH;
    FieldR3 savedE, savedH;
    FieldR3 nE, nH;
    // Jumps and fluxes: dim = (4*nfp, nK)
    FieldR3 dE, dH;
    FieldR3 dresE, dresH;
    // BC lists. nSMA, nPEC and nPMC are the number of BC of each kind.
    // BC are stored as pointers to memory positions in the jumps.
    // dim = (nK)
    UInt nSMA, nPEC, nPMC;
    UInt *SMAe, *SMAf, *PECe, *PECf, *PMCe, *PMCf;
    vector<DGSource*> source;
    vector<DGDispersive*> dispersive;
    void buildMaterials(
            const CellGroup& cells,
            const OptionsSolverDGTD& arg);
    void deduplicateVMaps(const CellGroup& cells);
    void allocateRHSAndJumps();
    void allocateMaps();
    void assignPointersToNeighbours(
            const CellGroup& cells,
            const Connectivities& map,
            const MeshVolume& mesh);
    void buildScalingFactors(const CellGroup& cells, const Connectivities& map);
    void buildEMSources(
            const EMSourceGroup& em,
            const BCGroup& bc,
            const Connectivities& maps,
            const CellGroup& cells);
    void BCToLocalArray(
            const BCGroup& bc,
            const CellGroup& cells,
            const Connectivities& map);
    vector<const BoundaryCondition*> removeNonLocalBCs(
            const CellGroup* cells,
            const vector<const BoundaryCondition*>& bc) const;
    bool checkPtrsToNeigh() const;
    void assignMatrices(const CellGroup& cells);
    void allocateFieldsForLTS();
    void buildCurvedFluxScalingFactors(const CellGroup& cells, const Connectivities& map);
};

inline UInt DGExplicit::getIndexOfElement(const UInt e) const {
    return (e * np);
}

#endif /* SOLVEREXPLICIT_H_ */
