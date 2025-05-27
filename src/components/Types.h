#pragma once

#include <array>
#include <vector>
#include <map>
#include <cfloat>
#include <mfem.hpp>

namespace maxwell {

const double TOLERANCE = 10.0 * std::numeric_limits<double>::epsilon();

struct FieldsForMovie {
	double Ex;
	double Ey;
	double Ez;
	double Hx;
	double Hy;
	double Hz;
};

using Time = double;
using FieldMovie = std::map<Time, double>;
using FieldMovies = std::map<Time,FieldsForMovie>;

using Point = std::vector<double>;
using Points = std::vector<Point>;

using FaceId = int;
using ElementId = int;
using BdrElementId = ElementId;
using El2Face = std::pair<ElementId, FaceId>;

using Attribute = int;
using BdrAttribute = Attribute;

enum FieldType { E, H };

enum SubMeshingMarkers {
	TotalFieldMarker = 1000,
	ScatteredFieldMarker = 2000,
	GlobalSubMeshMarker = 3000,
	NearToFarFieldMarker = 4000,
	SGBCMarker = 5000
};

enum class BdrCond {
	PEC,
	PMC,
	SMA,
	SurfaceCond,
	NearToFarField = 201,
	TotalFieldIn = 301,
	SGBC = 401
};

using InteriorFaceCoefficients = std::vector<double>;
using BdrFaceCoefficients = std::vector<double>;

using InteriorCoefficients = std::map<double, InteriorFaceCoefficients>;
using FluxBdrCoefficientsCentered = std::map<BdrCond, BdrFaceCoefficients>;
using FluxBdrCoefficientsUpwind = std::map<BdrCond, BdrFaceCoefficients>;
using FluxSrcCoefficientsCentered = std::map<BdrCond, BdrFaceCoefficients>;
using FluxSrcCoefficientsUpwind = std::map<BdrCond, BdrFaceCoefficients>;

using Direction = int;
static const Direction X{ 0 };
static const Direction Y{ 1 };
static const Direction Z{ 2 };

static const Attribute hesthavenMeshingTag{ 777 };


}