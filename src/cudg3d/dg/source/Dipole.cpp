#include "Dipole.h"
//
//DGDipole::DGDipole(
//      const Dipole& dip,
//      const BCGroup& bc,
//      const Connectivities& map,
//      const CellGroup& cells,
//      FieldR3& dE,
//      FieldR3& dH,
//      const Math::Int vmapM[faces][nfp])
//: Dipole(dip) {
//   initSource(bc, map, cells, dE, dH, vmapM);
//   // Determines total or scattered fields in the bc.
//   if (ETFNBe.size()) {
//      throw Error("Trying to set TF/SF in a not backed boundary.");
//   }
//   // Total field boundary.
//   vector<pair<size_t, size_t> > total;
//   total = getTotalFieldElemFaces(bc, map, cells);
//   tPos = new SphericalVector[ETInc.size() * nfp];
//   for (size_t i = 0; i < total.size(); i++) {
//      ElemId id = cells.getIdOfRelPos(total[i].first);
//      size_t f = total[i].second;
//      for (size_t j = 0; j < nfp; j++) {
//         tPos[i*nfp+j] =
//               cells.getPtrToCellWithId(id)->getSideNodePos(f,j) - position_;
//      }
//   }
//   // Scattered field boundary.
//   vector<pair<size_t,size_t> > scatt;
//   scatt = getScattFieldElemFaces(bc, map, cells);
//   sPos = new SphericalVector[ESInc.size() * nfp];
//   for (size_t i = 0; i < scatt.size(); i++) {
//      ElemId id = cells.getIdOfRelPos(scatt[i].first);
//      size_t f = scatt[i].second;
//      for (size_t j = 0; j < nfp; j++) {
//         sPos[i*nfp+j] =
//               cells.getPtrToCellWithId(id)->getSideNodePos(f,j) - position_;
//      }
//   }
//}
//
//void DGDipole::computeExcitation(
//      const Math::Real time,
//      const Math::Real minDT) {
//   computeExcitationField(ETInc, HTInc, tPos, ETInc.size(), time);
//   computeExcitationField(ESInc, HSInc, sPos, ESInc.size(), time);
//}
//
//void DGDipole::computeExcitationField(
//        FieldR3& EInc,
//        FieldR3& HInc,
//      const SphericalVector* vPos,
//      const size_t nE,
//      const Math::Real time) const {
//   // PURPOSE: Computes the dipole excitation.
//   // Chapter 2, R. Gomez's book. 2006.
//   // "Electromagnetic Field Theory for physicist and engineers.
//   // Section: Fields created by an infinitesimal current element with
//   // arbitrary time dependence.
//   // Uses the derivative of the gaussian.
//   // Otherwise charge accumulates and fields do not tend to zero.
//   Math::Real pos, pos2, pos3;
//   Math::Real expArg, expArg2;
//   Math::Real iT, iD;
//   Math::Real tDelayed, sint, cost;
//   SphericalVector sphE, sphH;
//   CVecR3 E, H;
//   // External field.
//   const size_t nFields = nfp * nE;
//   for (size_t j = 0; j < nFields; j++) {
//      pos = vPos[j].norm();
//      pos2 = pos * pos;
//      pos3 = pos2 * pos;
//      sint = sin(vPos[j].theta);
//      cost = cos(vPos[j].theta);
//      tDelayed = time - pos / Constants::c0; // Delayed time.
//      expArg = (tDelayed - gaussDelay_) / (spreadSqrt2_);
//      expArg2 = expArg * expArg;
//      iT = exp(- expArg2); // current @ delayed time.
//      iD = - iT * 2.0 * expArg / spreadSqrt2_; // derivative of current.
//      Math::Real iD2 = iD * (-2.0)* expArg / spreadSqrt2_
//            + iT * (-2.0) / spreadSqrt2_ / spreadSqrt2_;
//      Math::Real er = length_ * INV4PIEPS0 * 2.0 * cost
//            * (iT/pos3 + iD/(Constants::c0*pos2));
//      Math::Real et = length_ * INV4PIEPS0 * sint
//            * (iT/pos3
//                  + iD/(Constants::c0 * pos2)
//                  + iD2/(SPEED_OF_LIGHT_SQ*pos) );
//      Math::Real hp = length_ * INV4PI * sint
//            * (iD2/(pos*Constants::c0) + iD/pos2 );
//      // Spherical to Cartesian conversion.
//      E = vPos[j].convertSphericalVectorField(er, et, 0.0);
//      H = vPos[j].convertSphericalVectorField(0.0, 0.0, hp);
//      // use a rotation matrix applied on E and H.
//      EInc(x)[j] = E(0);
//      EInc(y)[j] = E(1);
//      EInc(z)[j] = E(2);
//      HInc(x)[j] = H(0);
//      HInc(y)[j] = H(1);
//      HInc(z)[j] = H(2);
//   }
//}
