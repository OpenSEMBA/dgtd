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
 * Arguments.cpp
 *
 *  Created on: Aug 24, 2012
 *      Author: luis
 */
#include "Options.h"

void OptionsSolverDGTD::initDefaults_() {
    timeIntegrator_ = lserk4;
    useMaxStageSizeForLTS_ = false;
    useLTS_ = true;
    growSmallerTiers_ = 0;
    maxNumberOfTiers_ = 0;
    upwinding_ = 1.0;
    PMLConstantConductivityProfile_ = false;
    PMLConductivity_ = 0.0;
}

OptionsSolverDGTD::OptionsSolverDGTD() {
    initDefaults_();
}

OptionsSolverDGTD::OptionsSolverDGTD(const OptionsSolver& base) :
        OptionsSolver(base) {
    initDefaults_();
}

OptionsSolverDGTD::~OptionsSolverDGTD() {

}

void OptionsSolverDGTD::set(const Arguments& arg) {
    OptionsSolver::set(arg);
    if (arg.has("timeintegrator")) {
        timeIntegrator_ = strToTimeIntegrator(arg.get("timeintegrator"));
    }
    if (arg.has("nolts")) {
        useLTS_ = false;
    }
    if (arg.has("usemaxstagesizeforlts")) {
        useMaxStageSizeForLTS_ = true;
    }
    if (arg.has("growsmallertiers")) {
        growSmallerTiers_ = atoi(arg.get("growsmallertiers").c_str());
    }
    if (arg.has("maxnumberoftiers")) {
        maxNumberOfTiers_ = atoi(arg.get("maxnumberoftiers").c_str());
    }
    if (arg.has("upwinding")) {
        upwinding_ = atof(arg.get("upwinding").c_str());
    }
    if (arg.has("pmluseconstantconductivity")) {
        PMLConstantConductivityProfile_ = true;
    }
    if (arg.has("pmlconductivity")) {
        PMLConductivity_ = atof(arg.get("pmlconductivity").c_str());
    }
}

void OptionsSolverDGTD::printInfo() const {
    OptionsSolver::printInfo();
    cout<< " -- Spatial discretization --" << endl;
    cout<< "Upwinding: " << upwinding_ << endl;
    cout<< " -- PML options -- " << endl;
    cout<< "Use Constant conductivity: "
            << PMLConstantConductivityProfile_ << endl;
    cout<< "Conductivity set to: " << PMLConductivity_ << endl;
    cout<< " -- Time integration --" << endl;
    cout<< "Integrator: ";
    switch (timeIntegrator_) {
    case lserk4:
        cout<< "4th Order Low-Storage Explicit Runge-Kutta" << endl;
        break;
    case verlet:
        cout<< "Verlet scheme" << endl;
        break;
    case lf2:
        cout<< "2nd Order Leapfrog (semi-defined)" << endl;
        break;
    case lf2full:
        cout<< "2nd Order Leapfrog (fully defined)" << endl;
        break;
    }
    cout<< "No LTS: " << useLTS_ << endl;
    cout<< "Use max. stage size: " << useMaxStageSizeForLTS_ << endl;
    cout<< "Grow smaller tiers: " << growSmallerTiers_ << endl;
    cout<< "Max. number of tiers: ";
    if (maxNumberOfTiers_ == 0) {
        cout<< "Unlimited" << endl;
    } else {
        cout<< maxNumberOfTiers_ << endl;
    }

}

OptionsSolverDGTD::TimeIntegrator OptionsSolverDGTD::getTimeIntegrator() const {
    return timeIntegrator_;
}

Real OptionsSolverDGTD::getUpwinding() const {
    return upwinding_;
}

void OptionsSolverDGTD::printHelp() const {
    cout<< " -i <name>" << endl;
    cout<< "     Specifies an input file."<< endl;
    cout<< " === General options ===" << endl;
    cout<< " -h or --help" << endl;
    cout<< "     Prints this help and exits." << endl;
    cout<< " === Time integrator options === " << endl;
    cout<< " --timeIntegrator <type> (defaults to lserk4)" << endl;
    cout<< "     lserk4  4th Order Low-Storage Explicit Runge-Kutta." << endl;
    cout<< "     verlet  2nd Order Verlet scheme." << endl;
    cout<< "     lf2     2nd Order Leapfrog scheme." << endl;
    cout<< "     lf2full 2nd Order LF. Completely defined states" << endl;
    cout<< " --noLTS " << endl;
    cout<< "     Deactivates Local Time Stepping." << endl;
    cout<< " --useMaxStageSizeForLTS" << endl;
    cout<< "     LTS will be done using max stage size of lserk4 as "
            << "ratio between time tiers. Default is (1 / Nstages). " << endl;
    cout<< " --timeStepSize <value from 0.0 to 1.0> (def. to 1.0)" << endl;
    cout<< "     1.0     Max. warranted stability." << endl;
    cout<< "    >1.0     Unsafe speed-up." << endl;
    cout<< " --maxNumberOfTiers <num>, --maxNumOfTiers <num>" << endl;
    cout<< "     Sets a maximum number of tiers for LTS." << endl;
    cout<< " --growSmallerTiers <num>" << endl;
    cout<< "     Includes <num> neighbours in smaller tiers." << endl;
    cout<< " === Spatial discretization options ===" << endl;
    cout<< " --upwinding <value from 0.0 to 1.0> (defaults to 1.0)" << endl;
    cout<< "     1.0     Full flux upwinding (Riemannian flux)." << endl;
    cout<< "     0.0     Centred flux." << endl;
    cout<< "  (0.0,1.0)  Partially penalized flux" << endl;
}

OptionsSolverDGTD::TimeIntegrator OptionsSolverDGTD::strToTimeIntegrator(
        const string& str) {
    if (!str.compare("verlet")) {
        return verlet;
    } else if (!str.compare("lf2")) {
        return lf2;
    } else if (!str.compare("lf2full")) {
        return lf2full;
    } else {
        return lserk4;
    }
}

UInt OptionsSolverDGTD::getGrowSmallerTiers() const {
    return growSmallerTiers_;
}

void OptionsSolverDGTD::setGrowSmallerTiers(UInt growSmallerTiers) {
    growSmallerTiers_ = growSmallerTiers;
}

UInt OptionsSolverDGTD::getMaxNumberOfTiers() const {
    return maxNumberOfTiers_;
}

void OptionsSolverDGTD::setMaxNumberOfTiers(UInt maxNumberOfTiers) {
    maxNumberOfTiers_ = maxNumberOfTiers;
}

Real OptionsSolverDGTD::getPMLConductivity() const {
    return PMLConductivity_;
}

void OptionsSolverDGTD::setPMLConductivity(Real pmlConductivity) {
    PMLConductivity_ = pmlConductivity;
}

bool OptionsSolverDGTD::isPMLConstantConductivityProfile() const {
    return PMLConstantConductivityProfile_;
}

void OptionsSolverDGTD::setPMLConstantConductivityProfile(
        bool pmlConstantConductivityProfile) {
    PMLConstantConductivityProfile_ = pmlConstantConductivityProfile;
}

void OptionsSolverDGTD::setTimeIntegrator(TimeIntegrator timeIntegrator) {
    timeIntegrator_ = timeIntegrator;
}

void OptionsSolverDGTD::setUpwinding(Real upwinding) {
    upwinding_ = upwinding;
}

bool OptionsSolverDGTD::isUseLTS() const {
    return useLTS_;
}

void OptionsSolverDGTD::setUseLTS(bool useLts) {
    useLTS_ = useLts;
}

bool OptionsSolverDGTD::isUseMaxStageSizeForLTS() const {
    return useMaxStageSizeForLTS_;
}

void OptionsSolverDGTD::setUseMaxStageSizeForLTS(bool useMaxStageSizeForLts) {
    useMaxStageSizeForLTS_ = useMaxStageSizeForLts;
}
