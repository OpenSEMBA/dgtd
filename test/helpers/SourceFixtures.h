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
	Sources res;
	res.push_back(std::move(std::make_unique<GaussianInitialField>(ft, d, spread, coeff, center)));
	return res;

}

static Sources buildRightTravelingWaveInitialField(const mfem::Vector& center)
{
	Sources res;
	res.push_back(std::move(std::make_unique<GaussianInitialField>(E, Y, 2.0, 1.0, center)));
	res.push_back(std::move(std::make_unique<GaussianInitialField>(H, Z, 2.0, 1.0, center)));
	return res;
}

static Sources buildSinusoidal(
	const FieldType& ft = E,
	const Direction& d = X,
	const std::vector<std::size_t> modes = { {1,0,0} },
	const double coefficient = 1.0,
	const mfem::Vector& cnt = mfem::Vector({ 0.5 }))
{
	Sources res;
	res.push_back(std::move(std::make_unique<SinusoidalInitialField>(ft, d, modes, coefficient, cnt)));
	return res;
}

}
}
}

