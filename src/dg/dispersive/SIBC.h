#pragma once

#include "Dispersive.h"
//#include "physicalModel/PMSurfaceSIBC.h"

//class DGSIBC : public DGDispersive, public PMSurfaceSIBC {
//   friend class PMSurface;
//public:
//   DGSIBC(
//         const PMSurfaceSIBC& mat_,
//         Math::Int ***map_,
//         const Math::Int vmapM[faces][nfp],
//         Math::Real ***ExP_,
//         Math::Real ***EyP_,
//         Math::Real ***EzP_,
//         Math::Real ***HxP_,
//         Math::Real ***HyP_,
//         Math::Real ***HzP_);
//   virtual ~DGSIBC();
//   void computeRHSElectricPolarizationCurrents(
//         const FieldR3& E,
//         const size_t e1,
//         const size_t e2);
//   void computeRHSMagneticPolarizationCurrents(
//         const FieldR3& H,
//         const size_t e1,
//         const size_t e2);
//   void computeRHSElectric(
//         FieldR3& rhsE,
//         const FieldR3& E,
//         const size_t e1, const size_t e2) const;
//   void computeRHSMagnetic(
//         FieldR3& rhsH,
//         const FieldR3& H,
//         const size_t e1, const size_t e2) const;
//   void addJumps(
//         FieldR3& dE, FieldR3& dH,
//         FieldR3& E, FieldR3& H,
//         const size_t e1, const size_t e2);
//   void addRHSToRes(
//         const size_t e1,
//         const size_t e2,
//         const Math::Real rka,
//         const Math::Real dt);
//   void updateWithRes(
//         const size_t e1,
//         const size_t e2,
//         const Math::Real rkb);
//private:
//   Math::Int ***map;
//   Math::Real ***ExP, ***EyP, ***EzP, ***HxP, ***HyP, ***HzP;
//   Math::Int vmapM[faces][nfp];
//   size_t nP;
//   size_t nE0, nED;
//   size_t *elem0, *elemD;
//   size_t *face0, *faceD;
//   CVecR3 *n0, *nD;
//   CVecR3 **Q0, **rhsQ0, **resQ0;
//   CVecR3 **QD, **rhsQD, **resQD;
//   CVecR3 *E0, *ED;
//   void computePolarizationFields(
//         const Math::Real *Hx, const Math::Real *Hy, const Math::Real *Hz,
//         const size_t e1, const size_t e2);
//};
//
//#endif /* SOLVERSIBC_H_ */
