#pragma once

#include "PML.h"
////
////class DGPMLMultiaxial : public DGPML {
////public:
////    DGPMLMultiaxial(
////            const PMVolumePML& mat,
////            const bool useConductivity,
////            const Math::Real conductivity);
////    virtual ~DGPMLMultiaxial();
////    void addRHSToRes(
////            const size_t e1, const size_t e2,
////            const Math::Real rka, const Math::Real dt);
////    void updateWithRes(
////            const size_t e1,
////            const size_t e2,
////            const Math::Real rkb);
////    virtual void computeRHSElectric(
////            FieldR3& rhsE,
////            const FieldR3& E,
////            const size_t e1, const size_t e2) const = 0;
////    virtual void computeRHSMagnetic(
////            FieldR3& rhsH,
////            const FieldR3& H,
////            const size_t e1, const size_t e2) const = 0;
////    virtual void computeRHSElectricPolarizationCurrents(
////            const FieldR3& E,
////            const size_t e1, const size_t e2) = 0;
////    virtual void computeRHSMagneticPolarizationCurrents(
////            const FieldR3& H,
////            const size_t e1, const size_t e2) = 0;
////protected:
////    FieldR3 J, M, resJ, resM, rhsJ, rhsM;
////};
////
//
//DGPMLMultiaxial::DGPMLMultiaxial(
//    const PMVolumePML& mat,
//    const CellGroup& cells,
//    const bool useConductivity,
//    const Math::Real conductivity) : DGPML(mat, cells) {
//    //    J.set(dof, 0.0);
//    //    resJ.set(dof, 0.0);
//    //    rhsJ.set(dof, 0.0);
//    //    M.set(dof, 0.0);
//    //    resM.set(dof, 0.0);
//    //    rhsM.set(dof, 0.0);
//}
//
//void DGPMLMultiaxial::addRHSToRes(
//    const size_t e1, const size_t e2,
//    const Math::Real rka, const Math::Real dt) {
//    //    size_t i, e;
//    //    for (i = 0; i < dof; i++) {
//    //        e = elem[(i / np) % nElem];
//    //        if (e1 <= e && e < e2) {
//    //            resJ[i] *= rka;
//    //            resJ[i] += rhsJ[i] * dt;
//    //            resM[i] *= rka;
//    //            resM[i] += rhsM[i] * dt;
//    //        }
//    //    }
//}
//
//void DGPMLMultiaxial::updateWithRes(
//    const size_t e1,
//    const size_t e2,
//    const Math::Real rkb) {
//    //    size_t i, e;
//    //#ifdef SOLVER_USE_OPENMP
//    //#pragma omp parallel for private(i, e)
//    //#endif
//    //    for (i = 0; i < dof; i++) {
//    //        e = elem[(i / np) % nElem];
//    //        if (e1 <= e && e < e2) {
//    //            J[i] += resJ[i] * rkb;
//    //            M[i] += resM[i] * rkb;
//    //        }
//    //    }
//}
