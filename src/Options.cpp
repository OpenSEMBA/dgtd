#include "Options.h"

namespace SEMBA {
namespace dgtd {

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
