#pragma once

#include "mfem.hpp"

namespace maxwell {

using namespace mfem;

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


}