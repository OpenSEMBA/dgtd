#pragma once

#include <array>
#include <vector>
#include <map>

namespace maxwell {

using Time = double;
using FieldFrame = std::vector<double>;
using FieldMovie = std::map<Time, FieldFrame>;

using Point = std::vector<double>;
using Points = std::vector<Point>;

enum FieldType {
	E,
	H
};

enum class FluxType {
	Centered,
	Upwind
};

struct FluxCoefficient {
	double beta;
};

enum class BdrCond {
	PEC,
	PMC,
	SMA,
	TotalField
};

struct MaxwellEvolOptions {
	FluxType fluxType{ FluxType::Upwind };
};


using Direction = int;
static const Direction X{ 0 };
static const Direction Y{ 1 };
static const Direction Z{ 2 };

enum class DisForm {
	Weak,
	Strong
};

enum class InitialFieldType {
	Gaussian,
	PlanarSinusoidal,
	PlaneWave
};


}