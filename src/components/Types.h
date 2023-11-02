#pragma once

#include <array>
#include <vector>
#include <map>
#include <mfem.hpp>

namespace maxwell {

struct GridFuncForFP {
	mfem::GridFunction Ex;
	mfem::GridFunction Ey;
	mfem::GridFunction Ez;
	mfem::GridFunction Hx;
	mfem::GridFunction Hy;
	mfem::GridFunction Hz;
};

struct FieldsForFP {
	double Ex;
	double Ey;
	double Ez;
	double Hx;
	double Hy;
	double Hz;
};

using Time = double;
using FieldMovie = std::map<Time, double>;
using FieldMovies = std::map<Time,FieldsForFP>;

using Point = std::vector<double>;
using Points = std::vector<Point>;

enum FieldType { E, H };
enum class FluxType { Centered, Upwind };

enum SubMeshingMarkers {
	TotalField = 1000,
	ScatteredField = 2000,
	Global_SubMesh = 3000
};

enum class BdrCond {
	NONE,
	PEC,
	PMC,
	SMA,
	TotalFieldIn = 301
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