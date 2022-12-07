#pragma once

#include "Dispersive.h"
#include "physicalModel/volume/Dispersive.h"

namespace SEMBA::dgtd::dg::dispersive {

class Volume : public Dispersive {
public:
//    DGDispersiveVolumic(
//            const PMVolumeDispersive&,
//            const CellGroup& cells);
//    void computeRHSElectricPolarizationCurrents(
//            const FieldR3& E,
//            const size_t e1, const size_t e2);
//    void computeRHSMagneticPolarizationCurrents(
//            const FieldR3& H,
//            const size_t e1, const size_t e2);
//    void computeRHSElectric(
//            FieldR3& rhsE,
//            const FieldR3& E,
//            const size_t e1, const size_t e2) const;
//    void computeRHSMagnetic(
//            FieldR3& rhsE,
//            const FieldR3& E,
//            const size_t e1, const size_t e2) const;
//    void addRHSToRes(
//            const size_t e1, const size_t e2,
//            const Math::Real rkb, const Math::Real dt);
//    void updateWithRes(
//            const size_t e1,
//            const size_t e2,
//            const Math::Real rkb);
//    void addJumps(
//            FieldR3& dE, FieldR3& dH,
//            FieldR3& E, FieldR3& H,
//            const size_t e1, const size_t e2);
//private:
//    static const size_t N = ORDER_N;
//    static const size_t np = (N+1) * (N+2) * (N+3) / 6;
//    //	PMVolumeDispersive mat;
//    size_t nElem;
//    size_t *elem;
//    size_t dof, drudeDof;
//    // Polarization currents. Size nK x Np x nPoles.
//    // Data is stored in nK vectors of Np components for each pole.
//    // First nK x Np data correspond to the first pole, and so on.
//    FieldC3 P, J;
//    FieldC3 rhsP, rhsJ;
//    FieldC3 resP, resJ;
};

}