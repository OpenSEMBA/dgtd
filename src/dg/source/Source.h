#pragma once

#include "core/source/Source.h"
#include "communications/Comm.h"
#include "boundaryConditions/Group.h"

namespace SEMBA {
namespace dgtd {
namespace dg {
namespace source {

class Source {
public:
    typedef enum {
        totalField,
        scatteredField,
        totalFieldNotBacked
    } BackingType;

    virtual ~DGSource() = default;
//    void addJumps(
//            const size_t e1,
//            const size_t e2);
//    virtual void computeExcitation(
//            const Math::Real intTime,
//            const Math::Real minDT) = 0;
//    virtual void printInfo() const = 0;
//protected:
//    const static size_t np = (N+1) * (N+2) * (N+3) / 6;
//    const static size_t np2 = np * 2;
//    const static size_t npnfp = np * nfp;
//    const static size_t npnp = np * np;
//    const static size_t nfpfaces = nfp * faces;
//    // Excitation fields.
//    FieldR3 ETInc, ESInc, EIncNB;
//    FieldR3 HTInc, HSInc, HIncNB;
//
//    vector<size_t> ETFe, ESFe, ETFNBe;
//    // Excitation total field jumps pointers.
//    Math::Real **dExT, **dEyT, **dEzT;
//    Math::Real **dHxT, **dHyT, **dHzT;
//    // Excitation scattered field jumps pointers.
//    Math::Real **dExS, **dEyS, **dEzS;
//    Math::Real **dHxS, **dHyS, **dHzS;
//    // Excitation total field not backed jumps.
//    Math::Real **dExTNB, **dEyTNB, **dEzTNB;
//    Math::Real **dHxTNB, **dHyTNB, **dHzTNB;
//    void initSource(
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells,
//            FieldR3& dE,
//            FieldR3& dH,
//            const Math::Int vmapM[faces][nfp]);
//    CVecR3* initPositions(
//            const vector<pair<size_t, size_t> >& elemFace,
//            const CellGroup& cells) const;
//    vector<pair<size_t, size_t>> getTotalFieldElemFaces(
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells) const;
//    vector<pair<size_t, size_t>> getScattFieldElemFaces(
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells) const;
//    vector<pair<size_t, size_t>> getTotalNotBackedFieldElemFaces(
//            const BCGroup& bc,
//            const Connectivities& map,
//            const CellGroup& cells) const;
};


}
}
}
}
