//// OpenSEMBA
//// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
////                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
////                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
////                    Daniel Mateos Romero            (damarro@semba.guru)
////
//// This file is part of OpenSEMBA.
////
//// OpenSEMBA is free software: you can redistribute it and/or modify it under
//// the terms of the GNU Lesser General Public License as published by the Free
//// Software Foundation, either version 3 of the License, or (at your option)
//// any later version.
////
//// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
//// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//// details.
////
//// You should have received a copy of the GNU Lesser General Public License
//// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
///*
// * DGWaveportRectangular.cpp
// *
// *  Created on: Aug 26, 2013
// *      Author: luis
// */
//
//#include "DGWaveguidePortRectangular.h"
//
//DGWaveportRectangular::DGWaveportRectangular(
//        const PortWaveguide& wp,
//        const MapGroup& map,
//        FieldR3& dE, FieldR3& dH,
//        const Int vmapM[faces][nfp]) :
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
//    //   Real zMax, zMin, yMax, yMin;
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
//    kcm = sqrt(pow((Real) getMode().first * Constants::pi/width, 2)
//            + pow((Real) getMode().second * Constants::pi/height, 2));
//    intrinsicImpedance = sqrt(Constants::mu0 / Constants::eps0);
//    gammaMSum = 0.0;
//}
//
//DGWaveportRectangular::~DGWaveportRectangular() {
//
//}
//
//void DGWaveportRectangular::computeExcitation(
//        const Real time,
//        const Real minDT) {
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
//void DGWaveportRectangular::printInfo() const {
//    cout << "DGWaveportRectangular::printInfo" << endl;
//    cout << "TO BE DONE." << endl;
//    // TODO DGWaveportRectangular::printInfo stub.
//}
//
//void DGWaveportRectangular::computeExcitationField(
//        FieldR3& EInc,
//        FieldR3& HInc,
//        const CVecR3* pos,
//        const size_t nE,
//        const Real time,
//        const Real minDT) {
//    //   const Real kcmSq = kcm * kcm;
//    //   //	Real f = getGauss(time, amplitude,delay,spread);
//    //   Real gamma0f = getMagnitude()->evaluate(time) / Constants::c0;
//    //   //	Real gammaMf =
//    //   //	 getNumericalGammaMGauss(time,minDT, amplitude,delay,spread, kcm);
//    //   if (excitationMode == Waveport::TE) {
//    //      const size_t nFields = nfp * nE;
//    //      const Real mConst = Constants::pi * getMode().first / width;
//    //      const Real nConst = Constants::pi * getMode().second / height;
//    //      for (size_t i = 0; i < nFields; i++) {
//    //         const Real yD = pos[i](1);
//    //         const Real zD = pos[i](2);
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
