#pragma once

#include <array>
#include <vector>
#include <map>

namespace maxwell {

using Time = double;
using FieldMovie = std::map<Time, double>;

using Point = std::vector<double>;
using Points = std::vector<Point>;

enum FieldType { E, H };
enum class FluxType { Centered, Upwind };

enum class BdrCond {
	NONE,
	PEC,
	PMC,
	SMA,
	TotalFieldIn = 301,
	TotalFieldOut = 302,
	TotalFieldInBacked = 303
};

using InteriorFaceCoefficient = std::vector<double>;
using BdrFaceCoefficient = std::vector<double>;

using InteriorCoefficients = std::map<FluxType, InteriorFaceCoefficient>;
using FluxCoefficientsCentered = std::map<BdrCond, BdrFaceCoefficient>;
using FluxCoefficientsUpwind = std::map<BdrCond, BdrFaceCoefficient>;

struct TFSFOrientationCoefficient {
	double orient;
};


struct MaxwellEvolOptions {
	FluxType fluxType{ FluxType::Upwind };
	bool spectral{ false };
	bool eigenvals{ false };
	bool marketFile{ false };
	int powerMethod{ 0 };
};


using Direction = int;
static const Direction X{ 0 };
static const Direction Y{ 1 };
static const Direction Z{ 2 };


}