#pragma once

//#include "DGDispersive.h"
//#include "physicalModel/PMVolumePML.h"
//
//class DGPML : public DGDispersive {
//public:
//    DGPML(const PMVolumePML& mat);
//    virtual ~DGPML();
//    void addJumps(
//            FieldR3& dE, FieldR3& dH,
//            FieldR3& E, FieldR3& H,
//            const size_t e1, const size_t e2);
//protected:
//    size_t dof;
//    size_t nElem;
//    size_t *elem;
//    bool useConstantConductivity;
//    static constexpr Math::Real sigDefault = 10e9;
//    Math::Real sig;
//    static const size_t N = ORDER_N;
//    static const size_t np = (N+1) * (N+2) * (N+3) / 6;
//    Math::Real **sig1, **sig2, **sig3;
//    Math::Real **sig11, **sig22, **sig33;
//    Math::Real **sig12, **sig23, **sig31;
//private:
//    void initConductivityMatrices(
//            const PMVolumePML& mat,
//            const CellGroup& cells);
//    void initConductivity(
//            Math::Real **sigma,
//            const size_t ori,
//            const PMVolumePML& mat,
//            const CellGroup& cells);
//};