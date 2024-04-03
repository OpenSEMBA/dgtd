#pragma once

#include "components/Sources.h"
#include "math/Calculus.h"

namespace maxwell {
namespace fixtures {
namespace sources {

static Sources buildGaussianInitialField(
	const FieldType& ft = E,
	const double spread = 0.1,
	const mfem::Vector& center_ = mfem::Vector({ 0.5 }),
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const int dimension = 1)
{
	mfem::Vector gaussianCenter(dimension);
	gaussianCenter = 0.0;
	
	Sources res;
	res.add(std::make_unique<InitialField>(
		Gaussian{ spread, gaussianCenter, dimension }, ft, p, center_)
	);
	return res;
}

static Sources buildResonantModeInitialField(
	const FieldType& ft = E,
	const Source::Polarization& p = Source::Polarization({ 0.0,0.0,1.0 }),
	const std::vector<std::size_t>& modes = { 1 })
{
	Sources res;
	Source::Position center_((int) modes.size());
	center_ = 0.0;
	res.add(
		std::make_unique<InitialField>(
			SinusoidalMode{ modes }, ft, p, center_
		)
	);
	return res;
}

static Sources buildGaussianPlanewave(
	double spread,
	double delay,
	const Source::Polarization& pol,
	const Source::Propagation& dir,
	const FieldType ft = FieldType::E
)
{
	Gaussian mag{ spread, mfem::Vector({-delay}) };
	Sources res;
	res.add(std::make_unique<Planewave>(mag, pol, dir, ft));
	return res;
}

static Sources buildPlanewaveInitialField(
	const Function& mf,
	const Source::Position& center_,
	const Source::Polarization& polIn,
	const Source::Propagation& propagationDir)
{
	Sources res;
	res.add(
			std::make_unique<InitialField>(mf, E, polIn, center_)
	);
	res.add(
			std::make_unique<InitialField>(mf, H, crossProduct(propagationDir, polIn), center_)
	);
	return res;
}

static Sources buildInitialField(
	const Function& mf)
{
	Sources res;
	res.add(
		std::make_unique<InitialField>(mf, E, Source::Polarization({0.0, 0.0, 1.0}), mfem::Vector({0.0,0.0,0.0}))
	);

	return res;
}

}
}
}

