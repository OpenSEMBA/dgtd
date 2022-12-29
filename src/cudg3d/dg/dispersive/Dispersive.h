#pragma once

//#include <cmath>

//#include "math/Constants.h"
//#include "math/Field.h"
//#include "physicalModel/PhysicalModel.h"
//#include "solver/dgtd/core/CellGroup.h"
//
//class DGDispersive {
//public:
//    static const size_t N = ORDER_N;
//    static const size_t nfp = (N+1) * (N+2) / 2;
//    static const size_t np = (N+1) * (N+2) * (N+3) / 6;
//    static const size_t faces = 4;
//    virtual ~DGDispersive() = 0;
//    virtual void computeRHSElectric(
//            FieldR3& rhsE,
//            const FieldR3& E,
//            const size_t e1, const size_t e2) const = 0;
//    virtual void computeRHSMagnetic(
//            FieldR3& rhsH,
//            const FieldR3& H,
//            const size_t e1, const size_t e2) const= 0;
//    virtual void computeRHSElectricPolarizationCurrents(
//            const FieldR3& E,
//            const size_t e1, const size_t e2) = 0;
//    virtual void computeRHSMagneticPolarizationCurrents(
//            const FieldR3& H,
//            const size_t e1, const size_t e2) = 0;
//    virtual void addRHSToRes(
//            const size_t e1, const size_t e2,
//            const Math::Real rka, const Math::Real dt) = 0;
//    virtual void updateWithRes(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rkb) = 0;
//    virtual void addJumps(
//            FieldR3& dE, FieldR3& dH,
//            FieldR3& E, FieldR3& H,
//            const size_t e1, const size_t e2) = 0;
//};
