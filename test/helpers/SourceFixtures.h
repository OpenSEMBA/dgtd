#pragma once

#include "maxwell/Sources.h"

namespace maxwell {
namespace fixtures {
namespace sources {

static Sources buildGaussianInitialField(
	const FieldType& ft = E,
	const Direction& d = X,
	const double spread = 0.1,
	const double coeff = 1.0,
	const mfem::Vector& center = mfem::Vector({ 0.5 }))
{
	return { GaussianInitialField(ft, d, spread, coeff, center) };
}

static Sources buildRightTravelingWaveInitialField(const mfem::Vector& center)
{
	return {
		GaussianInitialField(E, Y, 2.0, 1.0, center),
		GaussianInitialField(H, Z, 2.0, 1.0, center)
	};
}

}
}
}

