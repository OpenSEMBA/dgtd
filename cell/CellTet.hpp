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
//
//#include "cell/CellTet.h"
//
//namespace SEMBA {
//namespace Cudg3d {
//namespace Cell {
//
//template <int TET_N>
//CellTet<TET_N>::CellTet() {
//    base = NULL;
//}
//
//template <int TET_N>
//CellTet<TET_N>::~CellTet() {
//
//}
//
//template <int TET_N>
//void CellTet<TET_N>::init(
//        const Tet* base_,
//        const PMGroup* pMGroup) {
//    base = base_;
//    if (pMGroup != NULL) {
//        material = pMGroup->getId(base->getMatId())->castTo<PMVolumeClassic>();
//    } else {
//        material = NULL;
//    }
//    buildNodes();
//}
//
//template <int TET_N>
//array<typename CellTet<TET_N>::MatNpNp,3> CellTet<TET_N>::getCMatrices() const {
//    MatR44 cJ[SimplexTet<1>::ncp];
//    base->getCubatureJacobian(cJ);
//    double cJDet[SimplexTet<1>::ncp];
//    base->getCubatureJacobianDeterminant(cJDet, cJ);
//    MatNpNp invM = getMassMatrix(cJDet).invert();
//    StaMatrix<double,4,3> cJHat[SimplexTet<1>::ncp];
//    base->getCubatureJacobianHat(cJHat, cJ, cJDet);
//    array<MatNpNp,3> res;
//    for (size_t x = 0; x < 3; x++) {
//        res[x] = getCMatrix(x, invM, cJHat, cJ);
//    }
//    return res;
//}
//
//template <int TET_N>
//StaMatrix<double,TET_NP,TET_NP> CellTet<TET_N>::getConductivityWithGeometricProfile(
//        const PMVolumePML& mat,
//        const size_t type,
//        const double maxSigma) const {
////    static const size_t ncp = SimplexTet<TET_N>::ncp;
////    StaMatrix<double,np,np> res;
////    CVecR3 cPos[ncp];
////    base->getCubaturePositions(cPos);
////    int firstOrientation = mat.getFirstOrientationIndex();
////    int i = ((type % 10) - 1 + firstOrientation);
////    int j = ((type / 10) - 1 + firstOrientation);
////    const CVecR3 posBorder = mat.getPMLBeginningPosition();
////    const CVecR3 width = mat.getWidth();
////    // Computes sigma(x) = [(x - x_0)/width]^2 * maxSigma.
////    double sigma[ncp];
////    i = i % 3;
////    for (size_t c = 0; c < ncp; c++) {
////        sigma[c] = maxSigma;
////        sigma[c] *= pow((cPos[c](i) - posBorder(i))/width(i), 2);
////    }
////    // Computes second sigma if necessary.
////    if (j - firstOrientation >= 0) {
////        j = j % 3;
////        for (size_t c = 0; c < ncp; c++) {
////            sigma[c] *= maxSigma;
////            sigma[c] *= pow((cPos[c](j) - posBorder(j))/width(j), 2);
////        }
////    }
////    // Computes Mass Matrix.
////    res = getMassMatrixIntegratedWithScalar(sigma);
////    StaMatrix<double,np,np> invMM = getMassMatrix().invert();
////    return (invMM * res);
//}
//
//template <int TET_N>
//void CellTet<TET_N>::getCurvedLIFTnormal(
//        StaMatrix<double,np,nfp> LIFTn[3],
//        StaMatrix<double,np,nfp> LIFTcn[3],
//        StaMatrix<double,np,nfp> LIFTrn[3],
//        const size_t face) const {
//    throw ErrorNotImplemented("Not implemented");
//}
//
//template <int TET_N>
//size_t CellTet<TET_N>::getNbfp() const {
//    return base->numberOfSideCoordinates();
//}
//
//template <int TET_N>
//bool CellTet<TET_N>::isCurved() const {
//    return base->isCurved();
//}
//
//template <int TET_N>
//bool CellTet<TET_N>::isCurvedFace(const size_t f) const {
//    return base->isCurvedFace(f);
//}
//
//template <int TET_N>
//double CellTet<TET_N>::getVolume() const {
//    return base->getVolume();
//}
//
//template <int TET_N>
//const Tet* CellTet<TET_N>::getBase() const {
//    return base;
//}
//
//template <int TET_N>
//StaMatrix<double,TET_NP,TET_NP> CellTet<TET_N>::getCMatrix(
//        const size_t x,
//        const StaMatrix<double,np,np>& invM,
//        const StaMatrix<double,4,3> cJHat[SimplexTet<1>::ncp],
//        const StaMatrix<double,4,4> cJ[SimplexTet<1>::ncp]) const {
//    StaMatrix<double,TET_NP,TET_NP> res;
//    // Computes preliminary C matrices.
//    for (size_t i = 0; i < this->tet.np; i++) {
//        for (size_t j = 0; j < this->tet.np; j++) {
//            for (size_t k = 0; k < faces; k++) {
//                for (size_t c = 0; c < SimplexTet<1>::ncp; c++) {
//                    res(i,j) += this->tet.cwada[c][k](i,j) * cJHat[c](k,x);
//                }
//            }
//        }
//    }
//    // Multiplies by invMM.
//    res = invM * res;
//    res *= double(1.0 / 6.0);
//    return res;
//}
//
//template <int TET_N>
//StaMatrix<double,TET_NP,TET_NP> CellTet<TET_N>::getMassMatrix(
//        const double cJDet[SimplexTet<1>::ncp]) const{
//    StaMatrix<double,TET_NP,TET_NP> res;
//    for (size_t c = 0; c < SimplexTet<1>::ncp; c++) {
//        res += tet.cwaa[c] * cJDet[c];
//    }
//    res *= double(1.0 / 6.0);
//    return res;
//}
//
//template <int TET_N>
//StaMatrix<double,TET_NP,TET_NP>
//CellTet<TET_N>::getMassMatrix() const {
//    static const size_t ncp = SimplexTet<1>::ncp;
//    MatR44 cJ[SimplexTet<1>::ncp];
//    base->getCubatureJacobian(cJ);
//    double cJDet[ncp];
//    base->getCubatureJacobianDeterminant(cJDet, cJ);
//    return getMassMatrix(cJDet);
//}
//
//template <int TET_N>
//StaMatrix<double,TET_NP,TET_NP>
//CellTet<TET_N>::getMassMatrixIntegratedWithScalar(
//        const double cScalar[SimplexTet<1>::ncp]) const {
//    static const size_t ncp = SimplexTet<1>::ncp;
//    MatR44 cJ[SimplexTet<1>::ncp];
//    base->getCubatureJacobian(cJ);
//    double cJDetByScalar[ncp];
//    base->getCubatureJacobianDeterminant(cJDetByScalar, cJ);
//    for (size_t c = 0; c < ncp; c++) {
//        cJDetByScalar[c] *= cScalar[c];
//    }
//    return getMassMatrix(cJDetByScalar);
//}
//
//template <int TET_N>
//double CellTet<TET_N>::getAreaOfFace(size_t face) const {
//    return base->getAreaOfFace(face);
//}
//
//template <int TET_N>
//size_t CellTet<TET_N>::getNodeVertex(const size_t i) const {
//    return tet.vertex(i);
//}
//
//template <int TET_N>
//bool CellTet<TET_N>::isLocalSide(
//        const size_t face,
//        const SurfR* surf) const {
//    return (base->isLocalFace(face, *surf));
//}
//
//template <int TET_N>
//bool CellTet<TET_N>::isLocalSide(
//        const SurfR* surf) const {
//    for (size_t f = 0; f < faces; f++) {
//        if (isLocalSide(f, surf)) {
//            return true;
//        }
//    }
//    return false;
//}
//
//template <int TET_N>
//CVecR3 CellTet<TET_N>::getSideNodePos(
//        const size_t f, const size_t i) const {
//    return n[tet.sideNode(f, i)];
//}
//
//template <int TET_N>
//CVecR3 CellTet<TET_N>::getSideNormal(
//        const size_t f) const {
//    return base->getSideNormal(f);
//}
//
//template <int TET_N>
//size_t CellTet<TET_N>::getSideNode(
//        const size_t f, const size_t i) const {
//    return tet.sideNode(f, i);
//}
//
//template <int TET_N>
//const CoordR3* CellTet<TET_N>::getSideBaseNode(
//        const size_t f, const size_t i) const {
//    return base->getSideV(f,i);
//}
//
//template <int TET_N>
//const CoordR3* CellTet<TET_N>::getSideVertexBaseNode(size_t f, size_t i) const {
//    return (base->getSideVertex(f, i));
//}
//
//template <int TET_N>
//void CellTet<TET_N>::buildNodes() {
//    // Evaluates Lagrange's functions in positions specified by the
//    // simplex coordinates of tet.
//    double lagrEv[tet.np][base->numberOfCoordinates()];
//    for (size_t j = 0; j < tet.np; j++) {
//        for (size_t i = 0; i < base->numberOfCoordinates(); i++) {
//            lagrEv[j][i]= base->getTet().getLagr(i).eval(tet.coordinate(j));
//        }
//    }
//    // Computes nodes.
//    for (size_t j = 0; j < tet.np; j++) {
//        for (size_t i = 0; i < base->numberOfCoordinates(); i++) {
//            this->n[j] += *(base->getV(i)) * lagrEv[j][i];
//        }
//    }
//}
//
//template <int TET_N>
//bool CellTet<TET_N>::isFaceContainedInPlane(
//        const size_t face,
//        const CartesianPlane plane) const {
//    return base->isFaceContainedInPlane(face, plane);
//}
//
//template <int TET_N>
//void CellTet<TET_N>::printInfo() const {
//    cout << " --- CellTet information --- " << endl;
//    this->base->printInfo();
//    cout << "Nodes: " << this->np << endl;
//    for (size_t i = 0; i < this->np; i++) {
//        this->n[i].printInfo();
//        cout << endl;
//    }
//    cout << "Side Nodes and Normal vectors:" << endl;
//    for (size_t f = 0; f < base->numberOfFaces(); f++) {
//        cout << "Face " << f << ": " << endl;
//        for (size_t i = 0; i < this->nfp; i++) {
//            cout << "Pos: ";
//            this->getSideNodePos(f,i).printInfo();
//            cout << "Normal: ";
//            this->getSideNormal(f).printInfo();
//            cout << endl;
//        }
//    }
//    for (size_t f = 0; f < base->numberOfFaces(); f++) {
//        cout << "Area of face " << f << ": "
//                << this->base->getAreaOfFace(f) << endl;
//    }
//    cout << " Volume: " << this->base->getVolume() << endl;
//    cout << " Neighbour nodes: (vmapP)" << endl;
//    for (size_t f = 0; f < this->faces; f++) {
//        cout << " Face #" << f << ": ";
//        for (size_t i = 0; i < this->nfp; i++) {
//            cout << this->vmapP[f][i] << " ";
//        }
//        cout << endl;
//    }
//    cout << " Vertex info: " << endl;
//    for (size_t i = 0; i < this->vertices; i++) {
//        this->getV(i)->printInfo();
//        cout << endl;
//    }
//    cout << " Side Vertex Info: " << endl;
//    for (size_t f = 0; f < this->faces; f++) {
//        cout << "Face #" << f << endl;
//        for (size_t i = 0; i < this->nfp; i++) {
//            this->getSideVertexBaseNode(f,i)->printInfo();
//            cout << endl;
//        }
//    }
//    cout << " --- End of CellTet4 information --- " << endl;
//}
//
//template <int TET_N>
//CellTet4<TET_N>::CellTet4() {
//
//}
//
//template <int TET_N>
//CellTet4<TET_N>::~CellTet4() {
//
//}
//
//template <int TET_N>
//CellTet4<TET_N>::CellTet4(
//        const Tet* base_,
//        const PMGroup& pMGroup) {
//    this->init(base_, &pMGroup);
//}
//
//template <int TET_N>
//CellTet10<TET_N>::CellTet10() {
//
//}
//
//template <int TET_N>
//CellTet10<TET_N>::~CellTet10() {
//
//}
//
//template <int TET_N>
//CellTet10<TET_N>::CellTet10(
//        const Tet* base_,
//        const PMGroup& pMGroup) {
//    this->init(base_, &pMGroup);
//}
//
//template <int TET_N>
//void CellTet10<TET_N>::getCurvedLIFTnormal(
//        StaMatrix<double,np,nfp> LIFTn[3],
//        StaMatrix<double,np,nfp> LIFTcn[3],
//        StaMatrix<double,np,nfp> LIFTrn[3],
//        const size_t f) const {
//    // Side mass matrix.
//    vector<const CoordR3*> coord;
//    for (size_t i = 0; i < this->base->numberOfSideCoordinates(); i++) {
//        coord.push_back(this->base->getSideV(f,i));
//    }
//    Tri6 face(ElemId(0), &coord[0]);
//    static const size_t ncp = SimplexTri<2>::ncp;
//    double csdf[ncp];
//    face.getCubatureDifferentials(csdf);
//    CVecR3 csdn[ncp], csdcn[ncp], csdrn[ncp];
//    face.getCubatureNormals(csdn);
//    for (size_t c = 0; c < ncp; c++) {
//        csdcn[c](0) = csdn[c](0) * csdn[c](1);
//        csdcn[c](1) = csdn[c](1) * csdn[c](2);
//        csdcn[c](2) = csdn[c](2) * csdn[c](0);
//        csdrn[c] = (double) 1.0;
//        csdrn[c] -= csdn[c] * csdn[c];
//    }
//    StaMatrix<double,nfp,nfp> sMMn[3], sMMcn[3], sMMrn[3];
//    static const SimplexTri<TET_N> tri;
//    for (size_t i = 0; i < 3; i++) {
//        sMMn[i].zeros();
//        sMMcn[i].zeros();
//        sMMrn[i].zeros();
//        for (size_t c = 0; c < ncp; c++) {
//            sMMn[i] += tri.cwaa[c] * (csdf[c] / 2.0) * csdn[c](i);
//            sMMcn[i] += tri.cwaa[c] * (csdf[c] / 2.0) * csdcn[c](i);
//            sMMrn[i] += tri.cwaa[c] * (csdf[c] / 2.0) * csdrn[c](i);
//        }
//    }
//    // Mass matrix.
//    StaMatrix<int,nfp,np> RInt;
//    StaMatrix<int,np,nfp> RtrInt;
//    StaMatrix<double,np,nfp> R;
//    RInt = this->tet.R[f];
//    RtrInt = RInt.transpose();
//    R.assign(RtrInt);
//    StaMatrix<double,np,nfp> MMR;
//    MMR = this->getMassMatrix().invert() * R;
//    // Final form.
//    for (size_t i = 0; i < 3; i++) {
//        LIFTn[i] = MMR * sMMn[i];
//        LIFTcn[i] = MMR * sMMcn[i];
//        LIFTrn[i] = MMR * sMMrn[i];
//    }
//}
//
//}
//}
//}
