#pragma once

#include <array>
#include <vector>
#include <map>

#include "mfem.hpp"

namespace maxwell {

using namespace mfem;

using Time = double;
using CVec3 = std::array<double, 3>;
using FieldFrame = std::vector<CVec3>;
using FieldMovie = std::map<Time, FieldFrame>;

using Position = mfem::Vector;

enum FieldType {
	E,
	H
};

enum class FluxType {
	Centered,
	Upwind
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