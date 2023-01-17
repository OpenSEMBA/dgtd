#pragma once

#include "maxwell/Sources.h"

namespace maxwell {
namespace fixtures {
namespace sources {

static Sources buildGaussianInitialField(
	const FieldType& ft = E,
	const Direction& d = X,
	const double spread = 0.1,
	const mfem::Vector& center = mfem::Vector({ 0.5 }))
{
	Sources res;
	res.push_back(
		std::move(std::make_unique<InitialField>(
			GaussianFunction{ 1, spread, 1.0, center }, ft, d)
		)
	);
	return res;

}

static Sources buildRightTravelingWaveInitialField(const mfem::Vector& center)
{
	Sources res;
	GaussianFunction gauss{ 1, 2.0, 1.0, center };
	res.push_back(std::move(std::make_unique<InitialField>(gauss, E, Y)));
	res.push_back(std::move(std::make_unique<InitialField>(gauss, H, Z)));
	return res;
}

}
}
}

