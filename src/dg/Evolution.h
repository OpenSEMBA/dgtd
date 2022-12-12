#pragma once 

#include "model/Model.h"
#include "source/Group.h"

#include "VolumeModel.h"
#include "CellTet.h"
#include "Field.h"

namespace SEMBA::dgtd::dg {

using Model = SEMBA::Model::UnstructuredModel;
using EMSourceGroup = SEMBA::SourceGroup;

class Evolution {
public:
    using FaceNodeIndices = std::array<std::size_t, CellTet<POLYNOMIAL_ORDER>::nfp>;
	struct Options {
		Math::Real upwinding{ 1.0 };
	};
    Evolution(const VolumeModel&, const EMSourceGroup&, const Options&);
//    size_t getFieldDOFs();
//    const FieldR3& getRHSElectric() const;
//    const FieldR3& getRHSMagnetic() const;
private:
//    void computeRHS(const size_t e1, const size_t e2, const Math::Real localtime, const Math::Real rkdt);
//    void computeRHSElectric(const size_t e1, const size_t e2, const Math::Real localtime, const Math::Real minDT);
//    void computeRHSMagnetic(const size_t e1, const size_t e2, const Math::Real localtime, const Math::Real minDT);
//    void computeCurlsInRHSElectric(const size_t e1, const size_t e2);
//    void computeCurlsInRHSMagnetic(const size_t e1, const size_t e2);
//    void computeJumps(const size_t e1, const size_t e2, const Math::Real localTime, const Math::Real minDT);
//    void copyJumpsToResidueJumps(const size_t e1, const size_t e2);
//    void addFluxesToRHSElectric(const size_t e1, const size_t e2);
//    void addFluxesToRHSMagnetic(const size_t e1, const size_t e2);
//    void addFluxesToRHSElectric(const size_t e1, const size_t e2, const bool useResForUpw);
//    void addFluxesToRHSMagnetic(const size_t e1, const size_t e2, const bool useResForUpw);
//    void addStraightFluxesToRHSElectric(const size_t e1, const size_t e2, const bool useResForUpw);
//    void addStraightFluxesToRHSMagnetic(const size_t e1, const size_t e2, const bool useResForUpw);
//    void computePolarizationCurrentsRHS(const size_t e1, const size_t e2);
//    void computePolarizationCurrentsRHSElectric(const size_t e1, const size_t e2);
//    void computePolarizationCurrentsRHSMagnetic(const size_t e1, const size_t e2);
//    void addRHSToFieldsElectric(const size_t e1, const size_t e2, const Math::Real rkdt);
//    void addRHSToFieldsMagnetic(const size_t e1, const size_t e2, const Math::Real rkdt);
//    size_t getIndexOfElement(const size_t e) const;
//    void addRHSToResidueElectric(const size_t e1, const size_t e2, const Math::Real rkdt);
//    void addRHSToResidueMagnetic(const size_t e1, const size_t e2, const Math::Real rkdt);
//    void addRHSToRes(const size_t e1, const size_t e2, const Math::Real rka, const Math::Real dt);
//    void updateFieldsWithRes(const size_t e1, const size_t e2, const Math::Real rkb);
private:
    // - Maps.
    std::array<FaceNodeIndices, 4> vmapM;
    std::array<FaceNodeIndices, 16> vmapP;
//    Math::Int ***map_;
//    // Pointers to neighbour fields. dim = (nK, 4).
    Math::Real ***ExP, ***EyP, ***EzP, ***HxP, ***HyP, ***HzP;
	// Pointers to C. dim = (nK)
//    const Math::Real **Cx, **Cy, **Cz; 
//    // Fields and residuals: dim = (np,nK)
    FieldR3 rhsE, rhsH;
    //FieldR3 savedResE, savedResH;
    //FieldR3 savedE, savedH;
    //FieldR3 nE, nH;
//    // Jumps and fluxes: dim = (4*nfp, nK)
    FieldR3 dE, dH;
    //FieldR3 dresE, dresH;
//    // BC lists. nSMA, nPEC and nPMC are the number of BC of each kind.
//    // BC are stored as pointers to memory positions in the jumps.
//    // dim = (nK)
//    size_t nSMA, nPEC, nPMC;
//    size_t *SMAe, *SMAf, *PECe, *PECf, *PMCe, *PMCf;
//    vector<DGSource*> source;
//    vector<DGDispersive*> dispersive;
//    void buildMaterials(
//            const CellGroup& cells,
//            const OptionsSolverDGTD& arg);
//    void deduplicateVMaps(const CellGroup& cells);
    void allocateRHSAndJumps();
    void allocateMaps();
//    void assignPointersToNeighbours(
//            const CellGroup& cells,
//            const Connectivities& map,
//            const Mesh::Volume& mesh);
//    void buildScalingFactors(const CellGroup& cells, const Connectivities& map);
//    void buildEMSources(
//            const EMSourceGroup& em,
//            const BCGroup& bc,
//            const Connectivities& maps,
//            const CellGroup& cells);
//    void BCToLocalArray(
//            const BCGroup& bc,
//            const CellGroup& cells,
//            const Connectivities& map);
//    vector<const BoundaryCondition*> removeNonLocalBCs(
//            const CellGroup* cells,
//            const vector<const BoundaryCondition*>& bc) const;
//    bool checkPtrsToNeigh() const;
//    void assignMatrices(const CellGroup& cells);
//    void allocateFieldsForLTS();
//    void buildCurvedFluxScalingFactors(const CellGroup& cells, const Connectivities& map);
};

}