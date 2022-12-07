#pragma once 

namespace SEMBA::dgtd::dg {

//class Explicit {
//public:
//    Explicit(
//            const Mesh::Volume& mesh,
//            const PMGroup& pMGroup,
//            const EMSourceGroup& emsources,
//            const OptionsSolverDGTD& options,
//            Comm* comm);
//    virtual ~DGExplicit();
//    size_t getFieldDOFs();
//    const FieldR3& getRHSElectric() const;
//    const FieldR3& getRHSMagnetic() const;
//    void printInfo() const;
//protected:
//    void computeRHS(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real localtime,
//            const Math::Real rkdt);
//    void computeRHSElectric(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real localtime,
//            const Math::Real minDT);
//    void computeRHSMagnetic(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real localtime,
//            const Math::Real minDT);
//    void computeCurlsInRHSElectric(const size_t e1, const size_t e2);
//    void computeCurlsInRHSMagnetic(const size_t e1, const size_t e2);
//    void computeJumps(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real localTime,
//            const Math::Real minDT);
//    void copyJumpsToResidueJumps(
//            const size_t e1,
//            const size_t e2);
//    void addFluxesToRHSElectric(
//            const size_t e1, const size_t e2);
//    void addFluxesToRHSMagnetic(
//            const size_t e1, const size_t e2);
//    void addFluxesToRHSElectric(
//            const size_t e1, const size_t e2, const bool useResForUpw);
//    void addFluxesToRHSMagnetic(
//            const size_t e1, const size_t e2, const bool useResForUpw);
//    void addStraightFluxesToRHSElectric(const size_t e1, const size_t e2,
//            const bool useResForUpw);
//    void addStraightFluxesToRHSMagnetic(const size_t e1, const size_t e2,
//            const bool useResForUpw);
//    void addCurvedFluxesToRHSElectric(const size_t e1, const size_t e2,
//            const bool useResForUpw);
//    void addCurvedFluxesToRHSMagnetic(const size_t e1, const size_t e2,
//            const bool useResForUpw);
//    void computePolarizationCurrentsRHS(
//            const size_t e1, const size_t e2);
//    void computePolarizationCurrentsRHSElectric(
//            const size_t e1, const size_t e2);
//    void computePolarizationCurrentsRHSMagnetic(
//            const size_t e1, const size_t e2);
//    void addRHSToFieldsElectric(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rkdt);
//    void addRHSToFieldsMagnetic(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rkdt);
//    size_t getIndexOfElement(const size_t e) const;
//    void addRHSToResidueElectric(const size_t e1, const size_t e2,
//            const Math::Real rkdt);
//    void addRHSToResidueMagnetic(const size_t e1, const size_t e2,
//            const Math::Real rkdt);
//    void addRHSToRes(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rka,
//            const Math::Real dt);
//    void updateFieldsWithRes(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rkb);
//    void LTSSaveFieldsAndResidues(
//            const size_t fKSave,
//            const size_t lKSave);
//    void LTSLoadFieldsAndResidues(
//            const size_t fKSave,
//            const size_t lKSave);
//private:
//    // - Maps.
//    Math::Int vmapM[faces][nfp];
//    Math::Int vmapP[16][nfp];
//    Math::Int ***map_;
//    // Pointers to neighbour fields. dim = (nK, 4).
//    Math::Real ***ExP, ***EyP, ***EzP, ***HxP, ***HyP, ***HzP;
//    // Curved faces stuff ---------------------------------------------
//    size_t nCurvedFaces;
//    DGCurvedFace *curveFace;
//    const Math::Real **Cx, **Cy, **Cz; // Pointers to C. dim = (nK)
//    // Fields and residuals: dim = (np,nK)
//    FieldR3 rhsE, rhsH;
//    FieldR3 savedResE, savedResH;
//    FieldR3 savedE, savedH;
//    FieldR3 nE, nH;
//    // Jumps and fluxes: dim = (4*nfp, nK)
//    FieldR3 dE, dH;
//    FieldR3 dresE, dresH;
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
//    void allocateRHSAndJumps();
//    void allocateMaps();
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
//};

}