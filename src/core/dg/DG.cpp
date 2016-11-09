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
//#include "../dg/DG.h"
//
//DG::DG() {
//    nK = 0;
//    buildLIFT();
//}
//
//DG::~DG() {
//}
//
//void DG::setFieldsToRandom() {
//    static const Real min = -1.0;
//    static const Real max = 1.0;
//    E.setToRandom(min, max);
//    H.setToRandom(min, max);
//}
//
//void DG::setFieldsToGaussian(
//        const CellGroup& cells,
//        const Real amplitude,
//        CVecR3& polarization,
//        const CVecR3& gaussCenter,
//        const Real gaussWidth) {
//    CVecR3 aux;
//    Real expArg;
//    polarization.normalize();
//    for (size_t e = 0; e < nK; e++) {
//        ElemId id = cells.getIdOfRelPos(e);
//        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
//        for (size_t i = 0; i < np; i++) {
//            aux = cell->n[i] - gaussCenter;
//            expArg = aux.norm() / gaussWidth;
//            E.set(e*np + i, polarization*amplitude*exp(- expArg * expArg));
//        }
//    }
//    H.setAll((Real) 0.0);
//}
//
//void DG::setFieldsToHarmonics(
//        const CellGroup& cells,
//        const CVecI3& harmonics,
//        CVecR3& polarization) {
//    Real amp;
//    CVecR3 pos;
//    polarization.normalize();
//    for (size_t e = 0; e < nK; e++) {
//        ElemId id = cells.getIdOfRelPos(e);
//        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
//        for (size_t i = 0; i < np; i++) {
//            pos = cell->n[i];
//            if(harmonics(1) == 0) {
//                amp = sin(pos(0) * harmonics(0) * Constants::pi);
//            } else {
//                amp = sin(pos(0)*harmonics(0)*Constants::pi)
//				                     * sin(pos(1)*harmonics(1)*Constants::pi);
//            }
//            E.set(e*np + i, polarization*amp);
//        }
//    }
//    H.setAll((Real) 0.0);
//}
//
//void DG::setFieldsAndTimeFromResumeFile() {
//    //   string resumeFileName = "simulation.resume";
//    //   ifstream f_in(resumeFileName.c_str());
//    //   // Checks file existence.
//    //   if (!f_in) {
//    //      cerr << endl << "ERROR @ Solver::setFieldsFromResumeFile: "
//    //            << "File " << resumeFileName
//    //            << " does not exist, imposible to resume." << endl;
//    //   }
//    //   // Reads result header.
//    //   string trash;
//    //   getline(f_in, trash);
//    //   {
//    //      // Reads electric resume fields.
//    //      VectorModuleResult electricField(nK * np);
//    //      electricField.readResult(f_in);
//    //      // Copies result electric fields into the fast solver field vectors.
//    //      for (size_t i = 0; i < nK * np; i++) {
//    //         E.set(x)[i] = electricField.values[0][i];
//    //         E.set(y)[i] = electricField.values[1][i];
//    //         E.set(z)[i] = electricField.values[2][i];
//    //      }
//    //   }
//    //   {
//    //      // Reads magnetic resume fields.
//    //      VectorModuleResult magneticField(nK * np);
//    //      magneticField.readResult(f_in);
//    //      // Copies result electric fields into the fast solver field vectors.
//    //      for (size_t i = 0; i < nK * np; i++) {
//    //         H.set(x)[i] = magneticField.values[0][i];
//    //         H.set(y)[i] = magneticField.values[1][i];
//    //         H.set(z)[i] = magneticField.values[2][i];
//    //      }
//    //   }
//    //   f_in.close();
//    throw ErrorNotImplemented("Not implemented");
//}
//
//void DG::buildFieldScalingFactors(
//        const CellGroup& cells) {
//    oneOverEps = new Real[nK];
//    oneOverMu = new Real[nK];
//    for (size_t e = 0; e < nK; e++) {
//        ElemId id = cells.getIdOfRelPos(e);
//        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
//        oneOverEps[e] = 1.0 / (cell->material->getPermittivity());
//        oneOverMu[e]  = 1.0 / (cell->material->getPermeability());
//    }
//}
//
//void DG::buildFluxScalingFactors(
//        const CellGroup& cells,
//        const Connectivities& map) {
//    nAdm.setSize(nK*4);
//    nImp.setSize(nK*4);
//    rnAdm.setSize(nK*4);
//    rnImp.setSize(nK*4);
//    cnAdm.setSize(nK*4);
//    cnImp.setSize(nK*4);
//    // Straight faces -------------------------------------------------
//    for (size_t e = 0; e < nK; e++) {
//        ElemId id = cells.getIdOfRelPos(e);
//        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
//        Real impM, admM, impP, admP, impAv, admAv;
//        CVecR3 n, rn, cn;
//        CVecR3 nAdmAux, rnAdmAux, cnAdmAux;
//        CVecR3 nImpAux, rnImpAux, cnImpAux;
//        for (size_t f = 0; f < faces; f++) {
//            size_t i = e * faces + f;
//            // Computes Scaling factor.
//            Real fSc = 0.5 * cell->getAreaOfFace(f) / cell->getVolume();
//            // Computes local impedance and admittance..
//            impM = cell->material->getImpedance();
//            admM = cell->material->getAdmitance();
//            // Computes contiguous element impedance and admittance.
//            Face cellFace(cell->getBase(), f);
//            ElemId neighId;
//            if (map.isDomainBoundary(cellFace)) {
//                neighId = cellFace.first->getId();
//            } else {
//                neighId = map.getNeighFace(cellFace).first->getId();
//            }
//            const CellTet<ORDER_N>* neigh = cells.getPtrToCellWithId(neighId);
//            impP = neigh->material->getImpedance();
//            admP = neigh->material->getAdmitance();
//            impAv = (impM + impP) * 0.5;
//            admAv = (admM + admP) * 0.5;
//            // ----- Computes vectors ---------------------------------
//            n = cell->getSideNormal(f);
//            // Computes: rn = 1 - n .^ 2.
//            rn = Real(1.0);
//            rn -= n * n; // This is NOT cartesian scalar prod.
//            // Computes: cn = [ ny*nz  nz*nx nx*ny ].
//            cn(0) = n(0) * n(1);
//            cn(1) = n(1) * n(2);
//            cn(2) = n(2) * n(0);
//            // --- Assigns flux scaling factor ---
//            nAdm.set(i, n  * fSc * (admP / admAv));
//            nImp.set(i, n  * fSc * (impP / impAv));
//            cnAdm.set(i, cn * fSc / admAv);
//            cnImp.set(i, cn * fSc / impAv);
//            rnAdm.set(i, rn * fSc / admAv);
//            rnImp.set(i, rn * fSc / impAv);
//        }
//    }
//}
//
//void DG::init(
//        const OptionsSolverDGTD& options,
//        const PMGroup& pm,
//        const CellGroup& cells,
//        Comm* comm_) {
//    comm = comm_;
//    upwinding = options.getUpwinding();
//    nK = cells.getLocalSize();
//    buildMaterials(cells, options);
//    buildCMatrices(cells);
//    allocateFieldsAndRes();
//    resE.setAll((Real) 0.0);
//    resH.setAll((Real) 0.0);
//}
//
//void DG::addFluxesToRHS(
//        const size_t e1,
//        const size_t e2,
//        const Real localTime,
//        const Real minDT) {
//    computeJumps(e1,e2, localTime,minDT);
//    addFluxesToRHSElectric(e1,e2);
//    addFluxesToRHSMagnetic(e1,e2);
//}
//
//void DG::buildCMatrices(const CellGroup& cells) {
//    size_t e;
//#ifdef SOLVER_DEDUPLICATE_OPERATORS
//    for (e = 0; e < nK; e++) {
//        ElemId id = cells.getIdOfRelPos(e);
//        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
//        array<StaMatrix<Real,np,np>,3> C = cell->getCMatrices();
//        for (size_t i = 0; i < 3; i++) {
//            CList.insert(C[i]);
//        }
//    }
//#	else
//    CList = new StaMatrix<Real,np,np>[3*nK];
//#pragma omp parallel for private(e)
//    for (e = 0; e < nK; e++) {
//        size_t id = cells.getIdOfRelPos(e);
//        const CellTet<ORDER_N>* cell = cells.getPtrToCellWithId(id);
//        StaMatrix<Real,np,np> C[3];
//        cell->getCMatrices(C);
//        for (size_t i = 0; i < 3; i++) {
//            size_t j = 3 * e + i;
//            CList[j] = C[i];
//        }
//    }
//#	endif
//    assignMatrices(cells);
//}
//
//void DG::computeCurlsInRHS(
//        const size_t e1,
//        const size_t e2) {
//    computeCurlsInRHSElectric(e1, e2);
//    computeCurlsInRHSMagnetic(e1, e2);
//}
//
//void DG::updateFieldsWithResBase(
//        const size_t e1,
//        const size_t e2,
//        const Real rkb) {
//    size_t i, j, e;
//#pragma omp parallel for private(i,j,e)
//    for (e = e1; e < e2; e++) {
//        i = e * np;
//        for (j = 0; j < np; j++) {
//            E.set(x)[i] += resE(x)[i] * rkb;
//            E.set(y)[i] += resE(y)[i] * rkb;
//            E.set(z)[i] += resE(z)[i] * rkb;
//            H.set(x)[i] += resH(x)[i] * rkb;
//            H.set(y)[i] += resH(y)[i] * rkb;
//            H.set(z)[i] += resH(z)[i] * rkb;
//            i++;
//        }
//    }
//}
//
//
//void DG::copyFieldsInResidues(
//        const size_t e1,
//        const size_t e2) {
//    const size_t init = e1 * np;
//    const size_t end = e2 * np;
//    resE.copy(init, end, E);
//    resH.copy(init, end, H);
//}
//
//void DG::swapResiduesAndFields(
//        const size_t e1,
//        const size_t e2) {
//    const size_t init = e1 * np;
//    const size_t end = e2 * np;
//    E.swap(resE, init, end);
//    H.swap(resH, init, end);
//}
//
//void DG::buildScalingFactors(
//        const CellGroup& cells,
//        const Connectivities& map) {
//    buildFieldScalingFactors(cells);
//    buildFluxScalingFactors(cells, map);
//}
//
//void DG::buildLIFT() {
//    // Copies LIFT matrices into the solver matrix format.
//    // These are used only for linear elements.
//    static const SimplexTet<ORDER_N> tet;
//    StaMatrix<Real, np, nfp * 4> tmpLIFT;
//    for (size_t i = 0; i < np; i++) {
//        for (size_t f = 0; f < faces; f++) {
//            for (size_t j = 0; j < nfp; j++) {
//                tmpLIFT(i, f * nfp + j) = tet.LIFT[f](i, j);
//            }
//        }
//    }
//    tmpLIFT.convertToArray(MATRICES_ROW_MAJOR, LIFT);
//}
//
//void DG::allocateFieldsAndRes() {
//    size_t dof = getFieldDOFs();
//    E.setSize(dof/3);
//    H.setSize(dof/3);
//    resE.setSize(dof/3);
//    resH.setSize(dof/3);
//}
//
//const FieldR3* DG::getElectric() const {
//    return &E;
//}
//
//const FieldR3* DG::getMagnetic() const {
//    return &H;
//}
//
//size_t DG::getGlobalFieldPosOfVertex(pair<const ElemR*, size_t> vertex) const {
//    size_t e = getGlobalRelPosOfId(vertex.first->getId());
//    static const SimplexTet<ORDER_N> tet;
//    return (e*tet.np + tet.vertex(vertex.second));
//}
//
//vector<size_t> DG::getGlobalFieldPosOfFace(Face bound) const {
//    const size_t e = getGlobalRelPosOfId(bound.first->getId());
//    const size_t f = bound.second;
//    static const SimplexTet<ORDER_N> tet;
//    vector<size_t> res(tet.nfp, 0);
//    for (size_t i = 0; i < tet.nfp; i++) {
//        res[i] = e * tet.np + tet.sideNode(f,i);
//    }
//    return res;
//}
//
//vector<size_t> DG::getGlobalFieldPosOfVolume(const ElemId volId) const {
//    const size_t e = getGlobalRelPosOfId(volId);
//    static const SimplexTet<ORDER_N> tet;
//    vector<size_t> res(tet.np, 0);
//    for (size_t i = 0; i < tet.np; i++) {
//        res[i] = e * tet.np + i;
//    }
//    return res;
//}
