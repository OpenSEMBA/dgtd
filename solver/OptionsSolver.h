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
 * GlobalProblemData.h
 *
 *  Created on: Aug 27, 2012
 *      Author: luis
 */

#ifndef GLOBALPROBLEMDATA_H_
#define GLOBALPROBLEMDATA_H_

#include <cmath>
#include <iostream>
#include <utility>

using namespace std;

#include "Options.h"
#include "math/CartesianVector.h"

class OptionsSolver : public Options {
public:
    enum class EndingCondition {
        finalTime, numberOfTimeSteps
    };
    enum class Solver {
        ugrfdtd, cudg3d, none
    };
    OptionsSolver();
    virtual ~OptionsSolver();

    DEFINE_CLONE(OptionsSolver);

    void set(const Arguments& args);

    void setFinalTime(Real finalTime);
    void setSamplingPeriod(Real samplingPeriod);
    void setCFL(double cfl);
    void setNumberOfTimeSteps(UInt numberOfTimeSteps);
    void setFlush(Real flush);
    void setForceRestarting(bool forceRestarting);
    void setTimeStep(Real timeStep);
    void setResumeSimulation(bool resumeSimulation);
    void setAdditionalArguments(const string& additionalArguments);
    void setSolver(Solver solver);

    Real getFinalTime() const;
    Real getSamplingPeriod() const;
    Real getCFL() const;
    Real getTimeStep() const;
    Real getFlush() const;
    bool isForceRestarting() const;
    UInt getNumberOfTimeSteps() const;
    Solver getSolver() const;
    const string& getAdditionalArguments() const;
    bool isResumeSimulation() const;

    void printInfo() const;
    void printHelp() const;

    static string toStr(const OptionsSolver::Solver& solver);
    virtual string toArgsStr() const;

private:
    // Global
    Solver solver_;

    EndingCondition endingCondition_;
    Real ending_;
    Real timeStep_;
    Real cfl_;

    Real samplingPeriod_;
    bool forceRestarting_;
    bool resumeSimulation_;
    bool dontRun_;
    Real flush_;
    string additionalArguments_;
};

#endif /* GLOBALPROBLEMDATA_H_ */
