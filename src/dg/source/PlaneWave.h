#pragma once

//#include "../../dg/sources/DGSource.h"
//#include "sources/PlaneWave.h"
//
//using namespace std;
//
//class DGPlaneWave : public DGSource, public PlaneWave {
//public:
//    DGPlaneWave(
//            const PlaneWave& pw,
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells,
//            const Comm* comm,
//            FieldR3& dE, FieldR3& dH,
//            const Math::Int vmapM[faces][nfp]);
//    virtual ~DGPlaneWave();
//    void computeExcitation(
//            const Math::Real intTime,
//            const Math::Real minDT);
//    void printInfo() const;
//private:
//    Math::Real *kNPosTF; // dim(nETSF*SOLVER_NFP)
//    Math::Real *kNPosSF; // dim(nETSF*SOLVER_NFP)
//    Math::Real *kNPosTFNB; // dim(nETFNB*SOLVER_NFP)
//    void computeExcitationField(
//            FieldR3& EInc,
//            FieldR3& HInc,
//            const Math::Real* vPos,
//            const size_t nE,
//            const Math::Real intTime);
//    void initWaveNumberPosition(
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells,
//            const Comm* comm,
//            const Math::Int vmapM[faces][nfp]);
//};
//