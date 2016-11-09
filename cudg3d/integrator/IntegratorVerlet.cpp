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
// * SolverLeapfrog.cpp
// *
// *  Created on: Feb 22, 2013
// *      Author: luis
// */
//
//#include "../../dgtd/integrator/IntegratorVerlet.h"
//
//IntegratorVerlet::IntegratorVerlet() {
//}
//
//IntegratorVerlet::IntegratorVerlet(
//        const Mesh::Volume& mesh,
//        const PMGroup& pmGroup,
//        const OptionsSolverDGTD* arg) {
//    cfl_ = arg->getTimeStep();
//    cfl_ *= 0.9;
//    if (arg->getUpwinding() > 0.0) {
//        cfl_ *= 0.5;
//    }
//    init(mesh, pmGroup, arg);
//}
//
//IntegratorVerlet::~IntegratorVerlet() {
//
//}
//
//size_t IntegratorVerlet::getNumOfIterationsPerBigTimeStep(
//        const size_t e) const {
//    size_t nTiers = getNTiers();
//    size_t nStages = getNStages();
//    size_t tier = timeTierList_(e,1);
//    size_t iter = (nTiers - tier) * nStages;
//    return iter;
//}
//
//size_t IntegratorVerlet::getNStages() const {
//    return nStages;
//}
//
//Real IntegratorVerlet::getMaxTimeRatio() const {
//    return Real (0.5);
//}
//
//void IntegratorVerlet::timeIntegrate(
//        const Real time) const {
//    assert(solver != NULL);
//    if (doLTS) {
//        LTSTimeIntegration(time,getMaxDT(),getNTiers()-1);
//    } else {
//        updateFieldsVerlet(0,solver->nK,time,getMaxDT());
//    }
//}
//
//void IntegratorVerlet::LTSTimeIntegration(
//        Real localTime,
//        Real localdt,
//        const size_t tier) const {
//    size_t fK = getRange(tier, 0).first;
//    size_t lK = getRange(tier, 1).second;
//    if (tier > 0) {
//        LTSTimeIntegration(localTime, localdt/2.0, tier-1);
//    }
//    updateFieldsVerlet(fK,lK,localTime,localdt);
//    localTime += localdt / 2.0;
//    if (tier > 0) {
//        LTSTimeIntegration(localTime, localdt/2.0, tier-1);
//    }
//}
//
//void IntegratorVerlet::updateFieldsVerlet(
//        const size_t e1,
//        const size_t e2,
//        const Real localTime,
//        const Real rkdt) const {
//    solver->computeCurlsInRHSMagnetic(e1,e2);
//    solver->computeJumps(e1,e2,localTime,mindt);
//    solver->addFluxesToRHSMagnetic(e1,e2);
//    solver->addRHSToFieldsMagnetic(e1,e2,rkdt*0.5);
//    // ---
//    solver->computeCurlsInRHSElectric(e1,e2);
//    solver->computeJumps(e1,e2,localTime,mindt);
//    solver->addFluxesToRHSElectric(e1,e2);
//    solver->addRHSToFieldsElectric(e1,e2,rkdt);
//    // ---
//    solver->computeCurlsInRHSMagnetic(e1,e2);
//    solver->computeJumps(e1,e2,localTime,mindt);
//    solver->addFluxesToRHSMagnetic(e1,e2);
//    solver->addRHSToFieldsMagnetic(e1,e2,rkdt*0.5);
//}
