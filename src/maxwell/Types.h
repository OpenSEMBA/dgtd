#pragma once

namespace maxwell {

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