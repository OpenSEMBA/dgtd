// OpenSEMBA
// Copyright (C) 2015 Salvador Gonzalez Garcia        (salva@ugr.es)
//                    Luis Manuel Diaz Angulo         (lmdiazangulo@semba.guru)
//                    Miguel David Ruiz-Cabello Nuñez (miguel@semba.guru)
//                    Daniel Mateos Romero            (damarro@semba.guru)
//
// This file is part of OpenSEMBA.
//
// OpenSEMBA is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// OpenSEMBA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with OpenSEMBA. If not, see <http://www.gnu.org/licenses/>.
/*
 * GlobalProblemData.cpp
 *
 *  Created on: Aug 27, 2012
 *      Author: luis
 */

#include "OptionsSolver.h"

OptionsSolver::OptionsSolver () {
    // Global
    solver_ = Solver::none;
    endingCondition_ = EndingCondition::finalTime;
    ending_ = 0.0;
    samplingPeriod_ = 0.0;
    timeStep_ = 0.0;
    cfl_ = 0.8;
    forceRestarting_ = false;
    resumeSimulation_ = false;
    flush_ = 0.0;
    dontRun_ = false;
}

OptionsSolver::~OptionsSolver() {
}

void OptionsSolver::set(const Arguments& arg) {
    if (arg.has("h") || arg.has("help")) {
        printHelp();
    }
}

void OptionsSolver::printHelp() const {
    // TODO OptionsSolver printHelp
}

Real OptionsSolver::getFinalTime() const {
    if (endingCondition_ == EndingCondition::finalTime) {
        return ending_;
    } else {
        return ending_ * timeStep_;
    }
}

void OptionsSolver::setFinalTime(Real finalTime) {
    endingCondition_ = EndingCondition::finalTime;
    ending_ = finalTime;
}

Real OptionsSolver::getSamplingPeriod() const {
    return samplingPeriod_;
}

void OptionsSolver::setSamplingPeriod(Real samplingPeriod) {
    samplingPeriod_ = samplingPeriod;
}

Real OptionsSolver::getTimeStep() const {
    return timeStep_;
}

void OptionsSolver::setTimeStep(Real timeStep) {
    timeStep_ = timeStep;
}

double OptionsSolver::getCFL() const {
    return cfl_;
}

void OptionsSolver::setCFL(double cfl) {
    cfl_ = cfl;
}

UInt OptionsSolver::getNumberOfTimeSteps() const {
    if (endingCondition_ == EndingCondition::finalTime) {
        if (timeStep_ != 0.0) {
            return ceil(ending_ / timeStep_);
        } else {
            return (UInt) 0;
        }
    } else {
        return (UInt) ending_;
    }
}

void OptionsSolver::setNumberOfTimeSteps(UInt numberOfTimeSteps) {
    endingCondition_ = EndingCondition::numberOfTimeSteps;
    ending_ = (Real) numberOfTimeSteps;
}

void OptionsSolver::printInfo() const {
    cout<< " --- Solver parameters --- " << endl;
    cout<< "Solver: " << toStr(solver_) << endl;
    cout<< "Final time: " << getFinalTime() << endl;
    cout<< "Default sampling period: " << samplingPeriod_ << endl;
    cout<< "Time step: " << timeStep_ << endl;
}

const string& OptionsSolver::getAdditionalArguments() const {
    return additionalArguments_;
}

void OptionsSolver::setAdditionalArguments(const string& additionalArguments) {
    additionalArguments_ = additionalArguments;
}

Real OptionsSolver::getFlush() const {
    return flush_;
}

void OptionsSolver::setFlush(Real flush) {
    flush_ = flush;
}

bool OptionsSolver::isForceRestarting() const {
    return forceRestarting_;
}

void OptionsSolver::setForceRestarting(bool forceRestarting) {
    forceRestarting_ = forceRestarting;
}

bool OptionsSolver::isResumeSimulation() const {
    return resumeSimulation_;
}

void OptionsSolver::setResumeSimulation(bool resumeSimulation) {
    resumeSimulation_ = resumeSimulation;
}

string OptionsSolver::toStr(const OptionsSolver::Solver& solver) {
    switch (solver) {
    case Solver::cudg3d:
        return string("cudg3d");
    case Solver::ugrfdtd:
        return string("ugrfdtd");
    case Solver::none:
    default:
        return string("none");
    }
}

string OptionsSolver::toArgsStr() const {
    OptionsSolver defaultOptions;
    stringstream ss;
    if (getCFL() != defaultOptions.getCFL()) {
        ss << " -cfl " << getCFL();
    }
    if (isResumeSimulation()) {
        if (getTimeStep() != 0.0) {
            ss << " -r " << (UInt) floor(getFinalTime() / getTimeStep());
        } else {
            ss << " -r " << (UInt) floor(getFinalTime());
        }
    }
    if (isForceRestarting()) {
        ss << " -s";
    }
    if (getFlush() != defaultOptions.getFlush()) {
        ss << " -flush " << getFlush();
    }
    return ss.str();
}

OptionsSolver::Solver OptionsSolver::getSolver() const {
    return solver_;
}

void OptionsSolver::setSolver(OptionsSolver::Solver solver) {
    solver_ = solver;
}
