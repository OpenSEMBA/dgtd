#pragma once

//#include "Source.h"
//#include "sources/Dipole.h"

//class DGDipole : public DGSource {
//public:
//    DGDipole(
//            const Dipole& dip,
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells,
//            FieldR3& dE,
//            FieldR3& dH,
//            const Math::Int vmapM[faces][nfp]);
//    virtual ~DGDipole() = default;
//    void computeExcitation(
//            const Math::Real intTime,
//            const Math::Real minDT);
//    CVecR3 getMagnitude(const Math::Real time) const;

//private:
//    SphericalVector *tPos, *sPos;
//    void computeExcitationField(
//            FieldR3& EInc,
//            FieldR3& HInc,
//            const SphericalVector* vPos,
//            const size_t nE,
//            const Math::Real time) const;
//};
