#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

namespace maxwell {
namespace physicalConstants {
	constexpr double speedOfLight = 1.0;
	const double true_vacuum_permeability = 4.0 * M_PI * 1e-7;
	constexpr double trueSpeedOfLight = 299792458.0;
	constexpr double trueVacuumPermittivity = 8.8541878188e-12;
	const double trueFreeSpaceImpedance = 4.0 * M_PI * 1e-7 * 299792458.0;
}
}