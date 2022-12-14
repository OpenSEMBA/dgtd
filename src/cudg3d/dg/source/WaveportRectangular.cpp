#include "WaveportRectangular.h"
//
//DGWaveportRectangular::DGWaveportRectangular(
//        const PortWaveguide& wp,
//        const MapGroup& map,
//        FieldR3& dE, FieldR3& dH,
//        const Math::Int vmapM[faces][nfp]) :
//        PortWaveguide(wp) {
//    //   initSource(map, cells, dE, dH, vmapM);
//    //   // Computes positions.
//    //   vector<pair<size_t, size_t> > total;
//    //   total = getElemFaces(map, cells, totalField);
//    //   posTF = initPositions(total, cells);
//    //   if (!checkNormalsAreEqual(total, cells)) {
//    //      cerr << endl << "Total Normals are different" << endl;
//    //   }
//    //   vector<pair<size_t,size_t> > scatt;
//    //   scatt = getElemFaces(map, cells, scatteredField);
//    //   posSF = initPositions(scatt, cells);
//    //   if (!checkNormalsAreEqual(scatt, cells)) {
//    //      cerr << endl << "Scatt Normals are different" << endl;
//    //   }
//    //   vector<pair<size_t, size_t> > totalNB;
//    //   totalNB = getElemFaces(map, cells, totalFieldNotBacked);
//    //   posTFNB = initPositions(totalNB, cells);
//    //   if (!checkNormalsAreEqual(totalNB, cells)) {
//    //      cerr << endl << "Total Not Backed Normals are different" << endl;
//    //   }
//    //   // Compute waveport size.
//    //   Math::Real zMax, zMin, yMax, yMin;
//    //   if (nETF != 0) {
//    //      zMax = posTF[0](2);
//    //      zMin = posTF[0](2);
//    //      yMax = posTF[0](1);
//    //      yMin = posTF[0](1);
//    //   } else {
//    //      zMax = posTFNB[0](2);
//    //      zMin = posTFNB[0](2);
//    //      yMax = posTFNB[0](1);
//    //      yMin = posTFNB[0](1);
//    //   }
//    //   for (size_t i = 0; i < (nETF*nfp); i++) {
//    //      if (posTF[i](2) > zMax) {
//    //         zMax = posTF[i](2);
//    //      }
//    //      if (posTF[i](2) < zMin) {
//    //         zMin = posTF[i](2);
//    //      }
//    //      if (posTF[i](1) > yMax) {
//    //         yMax = posTF[i](1);
//    //      }
//    //      if (posTF[i](1) < yMin) {
//    //         yMin = posTF[i](1);
//    //      }
//    //   }
//    //   for (size_t i = 0; i < (nETFNB*nfp); i++) {
//    //      if (posTFNB[i](2) > zMax) {
//    //         zMax = posTFNB[i](2);
//    //      }
//    //      if (posTFNB[i](2) < zMin) {
//    //         zMin = posTFNB[i](2);
//    //      }
//    //      if (posTFNB[i](1) > yMax) {
//    //         yMax = posTFNB[i](1);
//    //      }
//    //      if (posTFNB[i](1) < yMin) {
//    //         yMin = posTFNB[i](1);
//    //      }
//    //   }
//    //   if (getSymXY() == Waveport::none) {
//    //      width = zMax - zMin;
//    //   } else {
//    //      width = (zMax - zMin) * 2.0;
//    //   }
//    //   if (getSymZX() == Waveport::none) {
//    //      height = yMax - yMin;
//    //   } else {
//    //      height = (yMax - yMin) * 2.0;
//    //   }
//    //   // Displaces origin to center of waveguide.
//    //   if (getSymXY() != Waveport::none) {
//    //      for (size_t i = 0; i < nETF*nfp; i++) {
//    //         posTF[i](2) += width / 2.0;
//    //      }
//    //      for (size_t i = 0; i < nESF*nfp; i++) {
//    //         posSF[i](2) += width / 2.0;
//    //      }
//    //      for (size_t i = 0; i < nETFNB*nfp; i++) {
//    //         posTFNB[i](2) += width / 2.0;
//    //      }
//    //   }
//    //   if (getSymZX() != Waveport::none) {
//    //      for (size_t i = 0; i < nETF*nfp; i++) {
//    //         posTF[i](1) += height / 2.0;
//    //      }
//    //      for (size_t i = 0; i < nESF*nfp; i++) {
//    //         posSF[i](1) += height / 2.0;
//    //      }
//    //      for (size_t i = 0; i < nETFNB*nfp; i++) {
//    //         posTFNB[i](2) += height / 2.0;
//    //      }
//    //   }
//    assert(false);
//#  warning "Waveports are not being ctrted."
//    // Stores modes.
//    excitationMode = getExcitationMode();
//    if (excitationMode != PortWaveguide::TE) {
//        cerr << endl << "ERROR @ DGWaveportRectangular" << endl;
//        cerr << endl << "Non TE mode not supported yet." << endl;
//        assert(false);
//        exit(-1);
//    }
//    // Computes kcm.
//    kcm = sqrt(pow((Math::Real) getMode().first * Constants::pi/width, 2)
//            + pow((Math::Real) getMode().second * Constants::pi/height, 2));
//    intrinsicImpedance = sqrt(Constants::mu0 / Constants::eps0);
//    gammaMSum = 0.0;
//}
//
//void DGWaveportRectangular::computeExcitation(
//        const Math::Real time,
//        const Math::Real minDT) {
//    //   computeExcitationField(
//    //         ExTInc, EyTInc, EzTInc, HxTInc, HyTInc, HzTInc,
//    //         posTF, nETF, time, minDT);
//    //   computeExcitationField(
//    //         ExSInc, EySInc, EzSInc, HxSInc, HySInc, HzSInc,
//    //         posSF, nESF, time, minDT);
//    //   computeExcitationField(
//    //         ExIncNB, EyIncNB, EzIncNB, HxIncNB, HyIncNB, HzIncNB,
//    //         posTFNB, nETFNB, time, minDT);
//}
//
//void DGWaveportRectangular::computeExcitationField(
//        FieldR3& EInc,
//        FieldR3& HInc,
//        const CVecR3* pos,
//        const size_t nE,
//        const Math::Real time,
//        const Math::Real minDT) {
//    //   const Math::Real kcmSq = kcm * kcm;
//    //   //	Math::Real f = getGauss(time, amplitude,delay,spread);
//    //   Math::Real gamma0f = getMagnitude()->evaluate(time) / Constants::c0;
//    //   //	Math::Real gammaMf =
//    //   //	 getNumericalGammaMGauss(time,minDT, amplitude,delay,spread, kcm);
//    //   if (excitationMode == Waveport::TE) {
//    //      const size_t nFields = nfp * nE;
//    //      const Math::Real mConst = Constants::pi * getMode().first / width;
//    //      const Math::Real nConst = Constants::pi * getMode().second / height;
//    //      for (size_t i = 0; i < nFields; i++) {
//    //         const Math::Real yD = pos[i](1);
//    //         const Math::Real zD = pos[i](2);
//    //         // --- Electric field, TE ---
//    //         //			ExInc[i] = 0.0;
//    //         EyInc[i] = gamma0f * intrinsicImpedance / kcmSq
//    //               * mConst * sin(mConst * zD) * cos(nConst * yD);
//    //         EzInc[i] = gamma0f * intrinsicImpedance / kcmSq
//    //               * nConst * cos(mConst * zD) * sin(nConst * yD);
//    //         // --- Magnetic field, TE ---
//    //         //			HxInc[i] = f * cos(mConst * zD) * cos(nConst * yD);
//    //         //			HyInc[i] = gammaMf / kcmSq
//    //         //			 * nConst * cos(mConst * zD) * sin(nConst * yD);
//    //         //			HzInc[i] = gammaMf / kcmSq
//    //         //			 * mConst * sin(mConst * zD) * cos(nConst * yD);
//    //      }
//    //   } else {
//    //      cerr << endl << "ERROR @ DGWaveportRectangular." << endl;
//    //      cerr << endl << "TE modes only." << endl;
//    //      assert(false);
//    //      exit(-1);
//    //   }
//}
