#include "Planewave.h"
//
//DGPlaneWave::DGPlaneWave(
//        const PlaneWave& pw,
//        const BCGroup& bc,
//        const Connectivities& map,
//        const CellGroup& cells,
//        const Comm* comm,
//        FieldR3& dE, FieldR3& dH,
//        const Math::Int vmapM[faces][nfp]) :
//        PlaneWave(pw) {
//    	initSource(bc, map, cells, dE, dH, vmapM);
//    	initWaveNumberPosition(bc, map, cells, comm, vmapM);
//}
//
//void DGPlaneWave::computeExcitation(
//        const Math::Real intTime, const Math::Real minDT) {
//    computeExcitationField(ETInc, HTInc, kNPosTF, ETInc.size(), intTime);
//    computeExcitationField(ESInc, HSInc, kNPosSF, ESInc.size(), intTime);
//    computeExcitationField(EIncNB, HIncNB, kNPosTFNB, EIncNB.size(), intTime);
//}
//
//void DGPlaneWave::computeExcitationField(
//        FieldR3& EInc,
//        FieldR3& HInc,
//        const Math::Real* vPos,
//        const size_t nE,
//        const Math::Real intTime) {
//    // Computes the plane wve excitation corresponding to the face f.
//    // The face contins nfp nodes and it is computed for time intTime.
//    const size_t nFields = nfp*nE;
//    for (size_t j = 0; j < nFields; j++) {
//        Math::Real delayedTime = intTime - vPos[j];
//        if (delayedTime >= 0) {
//            pair<CVecR3,CVecR3> EHInc = getElectromagneticField(delayedTime);
//            EInc(x)[j] = EHInc.first(0);
//            EInc(y)[j] = EHInc.first(1);
//            EInc(z)[j] = EHInc.first(2);
//            HInc(x)[j] = EHInc.second(0);
//            HInc(y)[j] = EHInc.second(1);
//            HInc(z)[j] = EHInc.second(2);
//        } else {
//            EInc(x)[j] = (Math::Real) 0.0;
//            EInc(y)[j] = (Math::Real) 0.0;
//            EInc(z)[j] = (Math::Real) 0.0;
//            HInc(x)[j] = (Math::Real) 0.0;
//            HInc(y)[j] = (Math::Real) 0.0;
//            HInc(z)[j] = (Math::Real) 0.0;
//        }
//    }
//}
//
//void DGPlaneWave::initWaveNumberPosition(
//        const BCGroup& bc,
//        const Connectivities& map,
//        const CellGroup& cells,
//        const Comm* comm,
//        const Math::Int vmapM[faces][nfp]) {
//    // Generates kNPosTSF for PW excitation.
//    // Computes krmin position.
//    Math::Real krmin = 0.0;
//    bool krminSet = false;
//    CVecR3 nodePos;
//    // Total field.
//    vector<pair<size_t, size_t> > total;
//    total = getTotalFieldElemFaces(bc, map, cells);
//    kNPosTF = new Math::Real[ETInc.size() * nfp];
//    for (size_t j = 0; j < ETInc.size(); j++) {
//        ElemId id = cells.getIdOfRelPos(total[j].first);
//        size_t f = total[j].second;
//        size_t pos = j * nfp;
//        for (size_t k = 0; k < nfp; k++) {
//            const CellTet<N>* cell = cells.getPtrToCellWithId(id);
//            nodePos = cell->n[vmapM[f][k]];
//            kNPosTF[pos + k] = getWaveDirection().dot(nodePos) / Constants::c0;
//            // Stores minimum kNPos value.
//            if (!krminSet) {
//                krminSet = true;
//                krmin = kNPosTF[pos + k];
//            }
//            else if (kNPosTF[pos + k] < krmin)
//                krmin = kNPosTF[pos + k];
//        }
//    }
//    // Scattered field.
//    vector<pair<size_t, size_t> > scatt;
//    scatt = getScattFieldElemFaces(bc, map, cells);
//    kNPosSF = new Math::Real[ESInc.size() * nfp];
//    for (size_t j = 0; j < ESInc.size(); j++) {
//        ElemId id = cells.getIdOfRelPos(scatt[j].first);
//        size_t f = scatt[j].second;
//        size_t pos = j * nfp;
//        for (size_t k = 0; k < nfp; k++) {
//            const CellTet<N>* cell = cells.getPtrToCellWithId(id);
//            nodePos = cell->n[vmapM[f][k]];
//            kNPosSF[pos + k] = getWaveDirection().dot(nodePos) / Constants::c0;
//            // Stores minimum kNPos value.
//            if (!krminSet) {
//                krminSet = true;
//                krmin = kNPosSF[pos + k];
//            }
//            else if (kNPosSF[pos + k] < krmin) {
//                krmin = kNPosSF[pos + k];
//            }
//        }
//    }
//    // Total field not backed.
//    vector<pair<size_t, size_t> > totalNotBacked;
//    totalNotBacked = getTotalNotBackedFieldElemFaces(bc, map, cells);
//    kNPosTFNB = new Math::Real[EIncNB.size() * nfp];
//    for (size_t j = 0; j < EIncNB.size(); j++) {
//        ElemId id = cells.getIdOfRelPos(totalNotBacked[j].first);
//        size_t f = totalNotBacked[j].second;
//        size_t pos = j * nfp;
//        for (size_t k = 0; k < nfp; k++) {
//            const CellTet<N>* cell = cells.getPtrToCellWithId(id);
//            nodePos = cell->n[vmapM[f][k]];
//            kNPosTFNB[pos + k] = getWaveDirection().dot(nodePos) / Constants::c0;
//            // Stores minimum kNPos value.
//            if (!krminSet) {
//                krminSet = true;
//                krmin = kNPosTFNB[pos + k];
//            }
//            else if (kNPosTFNB[pos + k] < krmin) {
//                krmin = kNPosTFNB[pos + k];
//            }
//        }
//    }
//    // Syncs minimum.
//    krmin = comm->reduceToGlobalMinimum(krmin);
//    // Adds krmin to kNPosTSF and kNPosTFNB.
//    for (size_t j = 0; j < ETInc.size() * nfp; j++) {
//        kNPosTF[j] -= krmin;
//    }
//    for (size_t j = 0; j < ESInc.size() * nfp; j++) {
//        kNPosSF[j] -= krmin;
//    }
//    for (size_t j = 0; j < EIncNB.size() * nfp; j++) {
//        kNPosTFNB[j] -= krmin;
//    }
//}
