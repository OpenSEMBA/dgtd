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

namespace SEMBA {
namespace Cudg3d {

void Options::initDefaults_() {
    timeIntegrator_ = TimeIntegrator::lserk4;
    useMaxStageSizeForLTS_ = false;
    useLTS_ = true;
    growSmallerTiers_ = 0;
    maxNumberOfTiers_ = 0;
    upwinding_ = 1.0;
    PMLConstantConductivityProfile_ = false;
    PMLConductivity_ = 0.0;
}

Options::Options() {
    initDefaults_();
}

Options::Options(const SEMBA::Solver::Options& base) :
        SEMBA::Solver::Options(base) {
    initDefaults_();
}

Options::~Options() {

}

void Options::addArguments(SEMBA::Argument::Group& arg) const {
    SEMBA::Solver::Options::addArguments(arg);
    arg.addOption(
            new Argument::Option<std::string>(
                    "Time integrator","timeintegrator")).choices(
                            {{"lserk4"}, {"lf2"}, {"lf2full"}, {"verlet"}});
    arg.addOption(new Argument::Switch("No LTS", "nolts"));
    arg.addOption(new Argument::Option<Math::Real>("Upwinding", "upwinding"));
}

void Options::set(const SEMBA::Solver::Settings& opt) {
    SEMBA::Solver::Options::set(opt);
    if (opt.existsName("Time integrator")) {
        setTimeIntegrator(
                strToTimeIntegrator(opt("Time integrator").getString()));
    }
    if (opt.existsName("No LTS")) {
        setUseLTS(!opt("No LTS").getBool());
    }
    if (opt.existsName("Upwinding")) {
        setUpwinding(opt("Upwinding").getReal());
    }
}

void Options::printInfo() const {
    SEMBA::Solver::Options::printInfo();
    cout << " -- Spatial discretization --" << endl;
    cout << "Upwinding: " << upwinding_ << endl;
    cout << "Time integration: " << toStr(getTimeIntegrator()) << endl;
    cout << "No LTS: " << useLTS_ << endl;
}

Options::TimeIntegrator Options::getTimeIntegrator() const {
    return timeIntegrator_;
}

Math::Real Options::getUpwinding() const {
    return upwinding_;
}

void Options::printHelp() const {
    SEMBA::Solver::Options::printHelp();
    cout<< " --timeIntegrator <type> (defaults to lserk4)" << endl;
    cout<< "     lserk4  4th Order Low-Storage Explicit Runge-Kutta." << endl;
    cout<< "     verlet  2nd Order Verlet scheme." << endl;
    cout<< "     lf2     2nd Order Leapfrog scheme." << endl;
    cout<< "     lf2full 2nd Order LF. Completely defined states" << endl;
    cout<< " --noLTS " << endl;
    cout<< "     Deactivates Local Time Stepping." << endl;
    cout<< " --upwinding <value from 0.0 to 1.0> (defaults to 1.0)" << endl;
    cout<< "     1.0     Full flux upwinding (Riemannian flux)." << endl;
    cout<< "     0.0     Centred flux." << endl;
    cout<< "  (0.0,1.0)  Partially penalized flux" << endl;
}

Options::TimeIntegrator Options::strToTimeIntegrator(const string& str) {
    if (!str.compare("verlet")) {
        return TimeIntegrator::verlet;
    } else if (!str.compare("lf2")) {
        return TimeIntegrator::lf2;
    } else if (!str.compare("lf2full")) {
        return TimeIntegrator::lf2full;
    } else {
        return TimeIntegrator::lserk4;
    }
}

size_t Options::getGrowSmallerTiers() const {
    return growSmallerTiers_;
}

void Options::setGrowSmallerTiers(size_t growSmallerTiers) {
    growSmallerTiers_ = growSmallerTiers;
}

size_t Options::getMaxNumberOfTiers() const {
    return maxNumberOfTiers_;
}

void Options::setMaxNumberOfTiers(size_t maxNumberOfTiers) {
    maxNumberOfTiers_ = maxNumberOfTiers;
}

Math::Real Options::getPMLConductivity() const {
    return PMLConductivity_;
}

void Options::setPMLConductivity(Math::Real pmlConductivity) {
    PMLConductivity_ = pmlConductivity;
}

bool Options::isPMLConstantConductivityProfile() const {
    return PMLConstantConductivityProfile_;
}

void Options::setPMLConstantConductivityProfile(
        bool pmlConstantConductivityProfile) {
    PMLConstantConductivityProfile_ = pmlConstantConductivityProfile;
}

void Options::setTimeIntegrator(TimeIntegrator timeIntegrator) {
    timeIntegrator_ = timeIntegrator;
}

void Options::setUpwinding(Math::Real upwinding) {
    upwinding_ = upwinding;
}

bool Options::isUseLTS() const {
    return useLTS_;
}

void Options::setUseLTS(bool useLts) {
    useLTS_ = useLts;
}

bool Options::isUseMaxStageSizeForLTS() const {
    return useMaxStageSizeForLTS_;
}

void Options::setUseMaxStageSizeForLTS(bool useMaxStageSizeForLts) {
    useMaxStageSizeForLTS_ = useMaxStageSizeForLts;
}

string Options::toStr(const Options::TimeIntegrator& timeIntegrator) {
    switch (timeIntegrator) {
    case TimeIntegrator::lserk4:
        return "4th Order Low-Storage Explicit Runge-Kutta";
    case TimeIntegrator::verlet:
        return "Verlet scheme";
    case TimeIntegrator::lf2:
        return "2nd Order Leapfrog (semi-defined)";
    case TimeIntegrator::lf2full:
        return "2nd Order Leapfrog (fully defined)";
    }

}

}
}
