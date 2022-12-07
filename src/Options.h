#pragma once 

#include "argument/Argument.h"
#include "solver/Options.h"

namespace SEMBA {
namespace dgtd {

class Options {
public:
	enum class TimeIntegrator {
		lserk4, verlet, lf2, lf2full
	};

    Math::Real upwinding_ = 1.0;
    TimeIntegrator timeIntegrator_ = TimeIntegrator::lserk4;
    bool useLTS_ = true;
    std::size_t growSmallerTiers_ = 0;
    std::size_t maxNumberOfTiers_ = 0;
    bool useMaxStageSizeForLTS_ = false;
    bool PMLConstantConductivityProfile_ = false;
    Math::Real PMLConductivity_ = 0.0;

	virtual ~Options() = default;

    void printHelp() const;

private:
	static TimeIntegrator strToTimeIntegrator(const string& str);
	static string toStr(const TimeIntegrator&);
};

}
}
