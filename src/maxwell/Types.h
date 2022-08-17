#pragma once

#include <array>
#include <vector>
#include <map>

namespace maxwell {

using Time = double;
using CVec3 = std::array<double, 3>;
using FieldFrame = std::vector<CVec3>;
using FieldMovie = std::map<Time, FieldFrame>;

enum FieldType {
	E,
	H
};

enum class FluxType {
	Centered,
	Upwind
};

struct FluxCoefficient {
	double alpha;
	double beta;
};

enum class BdrCond {
	PEC,
	PMC,
	SMA
};

enum Direction {
	X,
	Y,
	Z
};

enum class DisForm {
	Weak,
	Strong
};


}