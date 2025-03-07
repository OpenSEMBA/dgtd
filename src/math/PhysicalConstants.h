#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

namespace maxwell {
namespace physicalConstants {
	constexpr double speedOfLight = 1.0;
	constexpr double vacuumPermittivity = 1.0;
	constexpr double vacuumPermeability = 1.0;
	constexpr double vacuumPermeability_SI = 4.0 * M_PI * 1e-7;
	constexpr double speedOfLight_SI = 299792458.0;
	constexpr double vacuumPermittivity_SI = 8.8541878188e-12;
	constexpr double freeSpaceImpedance_SI = 4.0 * M_PI * 1e-7 * 299792458.0;
	constexpr double invFourPi = 1.0 / (4.0 * M_PI);
	constexpr double invFourPiEps0_SI = invFourPi / vacuumPermittivity_SI;
}
}