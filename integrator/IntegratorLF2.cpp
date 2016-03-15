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
//#include "IntegratorLF2.h"
//
//IntegratorLF2::IntegratorLF2() {
//}
//
//IntegratorLF2::IntegratorLF2(
// const Mesh::Volume& mesh,
// const PMGroup& pmGroup,
// const OptionsSolverDGTD* arg) {
//	cfl_ = arg->getTimeStep();
//	cfl_ *= 0.75;
//	if (arg->getUpwinding() > 0.0) {
//		cfl_ *= 0.5;
//	}
//	if (arg->isUseLTS()) {
//		cfl_ *= 0.95;
//	}
//	init(mesh, pmGroup, arg);
//}
//
//IntegratorLF2::~IntegratorLF2() {
//
//}
//
//void IntegratorLF2::setSolver(DG* solver_) {
//	solver = solver_;
//}
//
//void IntegratorLF2::timeIntegrate(
// const Real time) const {
//	assert(solver!=NULL);
//	Real dt = getMaxDT();
//	if (doLTS) {
//		LTSupdateFieldsElectric(time,dt,getNTiers()-1);
//		LTSupdateFieldsMagnetic(time,dt,getNTiers()-1);
//	} else {
//		updateFields(0,solver->nK,time,dt);
//	}
//}
//
//size_t IntegratorLF2::getNStages() const {
//	return nStages;
//}
//
//Real IntegratorLF2::getMaxTimeRatio() const {
//	return Real (1.0 / 3.0);
//}
//
//size_t
//IntegratorLF2::getNumOfIterationsPerBigTimeStep(
// const size_t e) const {
//	size_t nTiers = getNTiers();
//	size_t nStages = getNStages();
//	size_t tier = timeTierList_(e,1);
//	size_t iter = (nTiers - tier) * nStages;
//	return iter;
//}
//
//void IntegratorLF2::LTSupdateFieldsElectric(
// Real lTime,
// Real ldt,
// const size_t tier) const {
//	size_t e1 = getRange(tier, 0).first;
//	size_t e2 = getRange(tier, 1).second;
//	if (tier > 0) {
//		LTSupdateFieldsElectric(lTime, ldt/3.0, tier-1);
//		LTSupdateFieldsMagnetic(lTime+ldt/3.0,ldt/3.0, tier-1);
//		LTSupdateFieldsElectric(lTime+ldt*2.0/3.0, ldt/3.0, tier-1);
//	}
//	solver->computeRHSElectric(e1, e2, lTime, mindt);
//	solver->addRHSToFieldsElectric(e1, e2, ldt);
//}
//
//void IntegratorLF2::LTSupdateFieldsMagnetic(
// Real lTime,
// Real ldt,
// const size_t tier) const {
//	size_t fK = getRange(tier, 0).first;
//	size_t lK = getRange(tier, 1).second;
//	if (tier > 0) {
//		LTSupdateFieldsMagnetic(lTime, ldt/3.0, tier-1);
//		LTSupdateFieldsElectric(lTime+ldt/3.0,ldt/3.0, tier-1);
//		LTSupdateFieldsMagnetic(lTime+ldt*2.0/3.0, ldt/3.0, tier-1);
//	}
//	solver->computeRHSMagnetic(fK, lK, lTime, mindt);
//	solver->addRHSToFieldsMagnetic(fK, lK, ldt);
//}
//
//void IntegratorLF2::updateFields(
// const size_t e1,
// const size_t e2,
// const Real localTime,
// const Real rkdt) const {
//	solver->computeRHSElectric(e1, e2, localTime, mindt);
//	solver->addRHSToFieldsElectric(e1,e2,rkdt);
//	solver->computeRHSMagnetic(e1, e2, localTime, mindt);
//	solver->addRHSToFieldsMagnetic(e1,e2,rkdt);
//}
