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

using InteriorFaceCoefficients = std::vector<double>;
using BdrFaceCoefficients = std::vector<double>;

using InteriorCoefficients = std::map<FluxType, InteriorFaceCoefficients>;
using FluxBdrCoefficientsCentered = std::map<BdrCond, BdrFaceCoefficients>;
using FluxBdrCoefficientsUpwind = std::map<BdrCond, BdrFaceCoefficients>;
using FluxSrcCoefficientsCentered = std::map<BdrCond, BdrFaceCoefficients>;
using FluxSrcCoefficientsUpwind = std::map<BdrCond, BdrFaceCoefficients>;

struct TFSFOrientationCoefficient {
	double orient;
};

using Direction = int;
static const Direction X{ 0 };
static const Direction Y{ 1 };
static const Direction Z{ 2 };


}