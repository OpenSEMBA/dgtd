#pragma once

#include <array>
#include <vector>
#include <map>

namespace maxwell {

using Time = double;
using FieldMovie = std::map<Time, double>;

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

struct TFSFOrientationCoefficient {
	double orient;
};

enum class BdrCond {
	NONE,
	PEC,
	PMC,
	SMA,
	TotalFieldIn = 301,
	TotalFieldOut = 302
};

struct MaxwellEvolOptions {
	FluxType fluxType{ FluxType::Upwind };
	bool spectral{ false };
	bool eigenvals{ false };
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