#pragma once

namespace maxwell {

using Direction = std::size_t;

enum FieldType {
	Electric,
	Magnetic
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


}