#pragma once

#include "maxwell/Sources.h"
#include "maxwell/Calculus.h"

namespace maxwell {
namespace fixtures {
namespace sources {

static Sources buildGaussianInitialField(
	const FieldType& ft = E,
	const double spread = 0.1,
	const mfem::Vector& center = mfem::Vector({ 0.5 }),
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const int dimension = 1,
	const Source::CartesianAngles angles = Source::CartesianAngles({ 0.0,0.0,0.0 }))
{
	auto initialField{ 
		std::make_unique<InitialField>(
			Gaussian{ spread, dimension }, ft, p, center, angles) 
	};
	
	Sources res;
	res.push_back(std::move(initialField));
	return res;
}

static Sources buildResonantModeInitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const std::vector<std::size_t>& modes = { 1 },
	const int dim = 1)
{
	Sources res;
	Source::Position center(dim);
	center = 0.0;
	res.push_back(
		std::move(std::make_unique<InitialField>(
			SinusoidalMode{ dim, modes }, ft, p, center)
		)
	);
	return res;
}

static Sources buildPlaneWave(
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const double spread = 0.1,
	const double delay = 0.0,
	const int dimension = 1)
{
	Sources res;
	res.push_back(
		std::move(std::make_unique<PlaneWave>(
			TimeGaussian{ spread, delay, dimension }, p)
		)
	);
	return res;
}

static Sources buildPlanewaveInitialField(
	const MathFunction& mf,
	const Source::Position& center,
	const Source::Polarization& ePol,
	const mfem::Vector& propagationDir)
{
	Sources res;

	Source::Polarization hPol{ crossProduct(propagationDir, ePol)};

	res.push_back(std::move(std::make_unique<InitialField>(mf, E, ePol, center)));
	res.push_back(std::move(std::make_unique<InitialField>(mf, H, hPol, center)));

	return res;
}

static Sources buildConstantInitialField(
	const MathFunction& mf)
{
	Sources res;

	res.push_back(std::move(std::make_unique<InitialField>(mf, E, Source::Polarization({1.0, 0.0, 0.0}), mfem::Vector({0.0,0.0,0.0}))));
	res.push_back(std::move(std::make_unique<InitialField>(mf, E, Source::Polarization({0.0, 1.0, 0.0}), mfem::Vector({0.0,0.0,0.0}))));
	res.push_back(std::move(std::make_unique<InitialField>(mf, E, Source::Polarization({0.0, 0.0, 1.0}), mfem::Vector({0.0,0.0,0.0}))));
	res.push_back(std::move(std::make_unique<InitialField>(mf, H, Source::Polarization({1.0, 0.0, 0.0}), mfem::Vector({0.0,0.0,0.0}))));
	res.push_back(std::move(std::make_unique<InitialField>(mf, H, Source::Polarization({0.0, 1.0, 0.0}), mfem::Vector({0.0,0.0,0.0}))));
	res.push_back(std::move(std::make_unique<InitialField>(mf, H, Source::Polarization({0.0, 0.0, 1.0}), mfem::Vector({0.0,0.0,0.0}))));

	return res;
}

}
}
}

