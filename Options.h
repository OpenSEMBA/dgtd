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
 * Arguments.h
 *
 *  Created on: Aug 24, 2012
 *      Author: luis
 */

#ifndef ARGUMENTSDGTD_H_
#define ARGUMENTSDGTD_H_

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <omp.h>

using namespace std;
using namespace SEMBA;

#include "argument/Argument.h"
#include "solver/Options.h"

namespace Cudg3d {
namespace DGTD {

class Options : public Solver::Options {
public:
	enum class TimeIntegrator {
		lserk4, verlet, lf2, lf2full
	};

	Options();
	Options(const SEMBA::Solver::Options& base);
	virtual ~Options();

    DEFINE_CLONE(Options);

    void addArguments(Argument::Group& args) const;
    void set(const SEMBA::Solver::Settings& opts);

    void setGrowSmallerTiers(size_t growSmallerTiers);
    void setMaxNumberOfTiers(size_t maxNumberOfTiers);
    void setPMLConductivity(Math::Real pmlConductivity);
    void setPMLConstantConductivityProfile(bool);
    void setTimeIntegrator(TimeIntegrator timeIntegrator);
    void setUpwinding(Math::Real upwinding);
    void setUseLTS(bool useLts);
    void setUseMaxStageSizeForLTS(bool useMaxStageSizeForLts);

	TimeIntegrator getTimeIntegrator() const;
	Math::Real getUpwinding() const;
    size_t getGrowSmallerTiers() const;
    size_t getMaxNumberOfTiers() const;
    Math::Real getPMLConductivity() const;
    bool isPMLConstantConductivityProfile() const;
    bool isUseLTS() const;
    bool isUseMaxStageSizeForLTS() const;

	static TimeIntegrator strToTimeIntegrator(const string& str);
	static string toStr(const TimeIntegrator&);
	void printInfo() const;

private:
	Math::Real upwinding_;
	TimeIntegrator timeIntegrator_;
	bool useLTS_;
	size_t growSmallerTiers_;
	size_t maxNumberOfTiers_;
	bool useMaxStageSizeForLTS_;
	bool PMLConstantConductivityProfile_;
	Math::Real PMLConductivity_;

	void printHelp() const;
	void initDefaults_();

};

}
}

#endif /* ARGUMENTS_H_ */
