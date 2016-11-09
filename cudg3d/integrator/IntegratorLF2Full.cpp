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
//#include "../../dgtd/integrator/IntegratorLF2Full.h"
//
//IntegratorLF2Full::IntegratorLF2Full() {
//}
//
//IntegratorLF2Full::IntegratorLF2Full(
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
//
//}
//
//IntegratorLF2Full::~IntegratorLF2Full() {
//
//}
//
//size_t
//IntegratorLF2Full::getNStages() const {
//	return nStages;
//}
//
//Real
//IntegratorLF2Full::getMaxTimeRatio() const {
//	return ((Real)1.0 / (Real) nStages);
//}
//
//size_t
//IntegratorLF2Full::getNumOfIterationsPerBigTimeStep(const size_t e) const {
//	size_t nTiers = getNTiers();
//	size_t nStages = getNStages();
//	size_t tier = timeTierList_(e,1);
//	size_t stage = timeTierList_(e,2);
//	size_t iter = 0;
//	if (tier == 0) {
//		iter = (nTiers - tier) * nStages;
//	} else {
//		iter = (nTiers - tier) * (nStages + nStages - (stage+1));
//	}
//	return iter;
//}
//
//void
//IntegratorLF2Full::timeIntegrate(
// const Real time) const {
//	assert(solver != NULL);
//	if (doLTS) {
//		LTSTimeIntegration(time, getMaxDT(), getNTiers() - 1, 0);
//		LTSTimeIntegration(time, getMaxDT(), getNTiers() - 1, 0);
//	} else {
//		updateFields(0,solver->nK,time,getMaxDT());
//		updateFields(0,solver->nK,time,getMaxDT());
//	}
//}
//
//void
//IntegratorLF2Full::LTSTimeIntegration(
// Real localTime,
// Real localdt,
// const size_t tier,
// const size_t stage) const {
//	size_t fK, lK;
//	static const bool useResForUpwinding = true;
//	const size_t nStages = getNStages();
//	// Determines range of cells belonging to this tier and stage.
//	if (tier == getNTiers()-1) {
//		fK = getRange(tier, 0).first;
//		lK = getRange(tier, nStages-1).second;
//	} else {
//		fK = getRange(tier, 0).first;
//		lK = getRange(tier+1, stage).first;
//	}
//	// ------- Updates RHS in current tier ----------------------------
//	solver->computeCurlsInRHSElectric(fK,lK);
//	solver->computeCurlsInRHSMagnetic(fK,lK);
//	solver->computeJumps(fK,lK,localTime,mindt);
//	solver->addFluxesToRHSElectric(fK,lK,useResForUpwinding);
//	solver->addFluxesToRHSMagnetic(fK,lK,useResForUpwinding);
//	// ------- Recursively calls this function ------------------------
//	if (tier > 0) {
//		Real lts = localdt/((Real) nStages);
//		size_t lSaved = getRange(tier, nStages-2).second;
//		solver->LTSSaveFieldsAndResidues(fK,lSaved);
//		for (size_t s = 0; s < nStages; s++) {
//			LTSTimeIntegration(localTime + ((Real) s)*lts,
//			 lts, tier-1, nStages-s-1);
//		}
//		solver->LTSLoadFieldsAndResidues(fK,lSaved);
//	}
//	// ------- Updates field and residue data in current tier ---------
//	solver->addRHSToResidueElectric(fK,lK,localdt);
//	solver->addRHSToResidueMagnetic(fK,lK,localdt);
//	solver->copyJumpsToResidueJumps(fK,lK);
//	solver->swapResiduesAndFields(fK,lK);
//}
//
//void
//IntegratorLF2Full::updateFields(
// const size_t e1,
// const size_t e2,
// const Real localTime,
// const Real rkdt) const {
//	static const bool useResForUpwinding = true;
//	solver->computeCurlsInRHSElectric(e1,e2);
//	solver->computeCurlsInRHSMagnetic(e1,e2);
//	solver->computeJumps(e1,e2,localTime,mindt);
//	solver->addFluxesToRHSElectric(e1,e2,useResForUpwinding);
//	solver->addFluxesToRHSMagnetic(e1,e2,useResForUpwinding);
//	solver->addRHSToResidueElectric(e1,e2,rkdt);
//	solver->addRHSToResidueMagnetic(e1,e2,rkdt);
//	solver->copyJumpsToResidueJumps(e1,e2);
//	solver->swapResiduesAndFields(e1,e2);
//}
