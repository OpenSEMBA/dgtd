#include "PML.h"
//
//DGPML::DGPML(PMVolumePML& mat, CellGroup& cells) {
//    // Initializes Element list and dof.
//    vector<ElemId> auxList;
//    auxList.reserve(cells.getLocalSize());
//    for (size_t e = 0; e < cells.getLocalSize(); e++) {
//        const CellTet<ORDER_N>* cell = cells(e);
//        ElemId id = cell->getId();
//        if (cell->material->getId() == mat.getId()) {
//            auxList.push_back(cells.getRelPosOfId(id));
//        }
//        nElem = auxList.size();
//        elem = new size_t[nElem];
//        for (size_t j = 0; j < nElem; j++) {
//            elem[j] = auxList[j];
//        }
//    }
//    dof = np * nElem;
//    // Initializes conductivity matrices.
//    if (!useConstantConductivity) {
//        initConductivityMatrices(mat, cells);
//    }
//}
//
//void DGPML::initConductivity(
//        Math::Real **sigma,
//        const size_t orientation,
//        const PMVolumePML& mat,
//        const CellGroup& cells) {
//    StaMatrix<Math::Real,np,np> auxSig;
//    for (size_t e = 0; e < nElem; e++) {
//        ElemId id = cells.getIdOfRelPos(elem[e]);
//        const CellTet<ORDER_N>* cell;
//        cell = cells.getPtrToCellWithId(id);
//        auxSig =
//                cell->getConductivityWithGeometricProfile(mat, orientation, sig);
//        auxSig.convertToArray(MATRICES_ROW_MAJOR, sigma[e]);
//    }
//}
//
//void DGPML::initConductivityMatrices(
//        const PMVolumePML& mat,
//        const CellGroup& cells) {
//    assert(elem != NULL);
//    if (mat.isUniaxial()) {
//        sig1 = new Math::Real*[nElem];
//        sig11 = new Math::Real*[nElem];
//        for (size_t e = 0; e < nElem; e++) {
//            sig1[e] = new Math::Real[np*np];
//            sig11[e] = new Math::Real[np*np];
//        }
//        initConductivity(sig1, 1, mat, cells);
//        initConductivity(sig11, 11, mat, cells);
//    } else if (mat.isBiaxial()) {
//        sig1 = new Math::Real*[nElem];
//        sig2 = new Math::Real*[nElem];
//        sig11 = new Math::Real*[nElem];
//        sig22 = new Math::Real*[nElem];
//        sig12= new Math::Real*[nElem];
//        for (size_t e = 0; e < nElem; e++) {
//            sig1[e] = new Math::Real[np*np];
//            sig2[e] = new Math::Real[np*np];
//            sig11[e] = new Math::Real[np*np];
//            sig22[e] = new Math::Real[np*np];
//            sig12[e] = new Math::Real[np*np];
//        }
//        initConductivity(sig1, 1, mat, cells);
//        initConductivity(sig2, 2, mat, cells);
//        initConductivity(sig11, 11, mat, cells);
//        initConductivity(sig22, 22, mat, cells);
//        initConductivity(sig12, 12, mat, cells);
//    } else {
//        assert(mat.isTriaxial());
//        sig1 = new Math::Real*[nElem];
//        sig2 = new Math::Real*[nElem];
//        sig3 = new Math::Real*[nElem];
//        sig11 = new Math::Real*[nElem];
//        sig22 = new Math::Real*[nElem];
//        sig33 = new Math::Real*[nElem];
//        sig12= new Math::Real*[nElem];
//        sig23= new Math::Real*[nElem];
//        sig31= new Math::Real*[nElem];
//        for (size_t e = 0; e < nElem; e++) {
//            sig1[e] = new Math::Real[np*np];
//            sig2[e] = new Math::Real[np*np];
//            sig3[e] = new Math::Real[np*np];
//            sig11[e] = new Math::Real[np*np];
//            sig22[e] = new Math::Real[np*np];
//            sig33[e] = new Math::Real[np*np];
//            sig12[e] = new Math::Real[np*np];
//            sig23[e] = new Math::Real[np*np];
//            sig31[e] = new Math::Real[np*np];
//        }
//        initConductivity(sig1, 1, mat, cells);
//        initConductivity(sig2, 2, mat, cells);
//        initConductivity(sig3, 3, mat, cells);
//        initConductivity(sig11, 11, mat, cells);
//        initConductivity(sig22, 22, mat, cells);
//        initConductivity(sig33, 33, mat, cells);
//        initConductivity(sig12, 12, mat, cells);
//        initConductivity(sig23, 23, mat, cells);
//        initConductivity(sig31, 31, mat, cells);
//    }
//}
//
//void DGPML::addJumps(
//        FieldR3& dE, FieldR3& dH,
//        FieldR3& E, FieldR3& H,
//        const size_t e1, const size_t e2) {
//}
