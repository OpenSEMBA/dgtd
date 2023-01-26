#pragma once

#include "maxwell/Sources.h"

namespace maxwell {
namespace fixtures {
namespace sources {

static Sources buildGaussianInitialField(
	const FieldType& ft = E,
	const Direction& d = X,
	const double spread = 0.1,
	const mfem::Vector& center = mfem::Vector({ 0.5 }),
	const int dimension = 1)
{
	Sources res;
	res.push_back(
		std::move(std::make_unique<InitialField>(
			Gaussian{ dimension, spread, 1.0, center }, ft, d)
		)
	);
	return res;

}

static Sources buildPlaneWave(
	const Direction& d = X,
	const double spread = 0.1,
	const double delay = 0.0,
	const double normalization = 1.0,
	const int dimension = 1)
{
	Sources res;
	res.push_back(
		std::move(std::make_unique<PlaneWave>(
			TimeGaussian{ dimension, spread, normalization, delay }, d)
		)
	);
	return res;
}

static Sources buildRightTravelingWaveInitialField(const Gaussian& gauss)
{
	Sources res;
	res.push_back(std::move(std::make_unique<InitialField>(gauss, E, Y)));
	res.push_back(std::move(std::make_unique<InitialField>(gauss, H, Z)));
	return res;
}

}
}
}

