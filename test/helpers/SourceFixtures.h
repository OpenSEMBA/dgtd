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
	const double spread = 0.1,
	const double delay = 0.0,
	const int dimension = 1,
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const FieldType ft = E,
	const mfem::Vector& propagationDir = mfem::Vector({1.0,0.0,0.0}),
	const Source::Position& center = Source::Position({0.0,0.0,0.0}),
	const Source::CartesianAngles& angles = Source::CartesianAngles({0.0,0.0,0.0}))
{
	auto altFt{ ft == E ? H : E };
	Sources res;
	res.push_back(
		std::move(std::make_unique<TimeVaryingField>(
			TimeGaussian{ spread, delay, dimension }, 
			p, 
			ft,
			center, 
			angles)
		)
	);
	res.push_back(
		std::move(std::make_unique<TimeVaryingField>(
			TimeGaussian{ spread, delay, dimension },
			crossProduct(propagationDir, p),
			altFt,
			center,
			angles)));
	return res;
}

static Sources buildPlanewaveInitialField(
	const MathFunction& mf,
	const FieldType& ft,
	const Source::Position& center,
	const Source::Polarization& polIn,
	const mfem::Vector& propagationDir,
	const Source::CartesianAngles& angles = Source::CartesianAngles({0.0,0.0,0.0}))
{
	Sources res;

	auto altFt{ ft == E ? H : E };

	res.push_back(std::move(std::make_unique<InitialField>(mf,    ft,                               polIn, center, angles)));
	res.push_back(std::move(std::make_unique<InitialField>(mf, altFt, crossProduct(propagationDir, polIn), center, angles)));

	return res;
}

static Sources buildConstantInitialField(
	const MathFunction& mf)
{
	Sources res;

	res.push_back(std::move(std::make_unique<InitialField>(mf, E, Source::Polarization({0.0, 0.0, 1.0}), mfem::Vector({0.0,0.0,0.0}))));

	return res;
}

}
}
}

