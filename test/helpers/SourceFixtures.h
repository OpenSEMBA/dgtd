#pragma once

#include "maxwell/Sources.h"

namespace maxwell {
namespace fixtures {
namespace sources {

	static Sources buildGaussianInitialField(
		const FieldType& ft = E,
		const double spread = 0.1,
		const mfem::Vector& center = mfem::Vector({ 0.5 }),
		const Source::Polarization& p = Source::Polarization{ 0.0,0.0,1.0 },
		const int dimension = 1,
		const double rotAngle = 0.0)
{
	assert(center.Size() >= dimension);
	mfem::Vector gaussianCenter(dimension);
	if (center.Size() == dimension) {
		gaussianCenter = center;
	}
	else
	{
		gaussianCenter[0] = center.Norml2();
	}
	
	auto initialField{ 
		std::make_unique<InitialField>(Gaussian{ spread, dimension }, ft, p, center, rotAngle) };
	
	Sources res;
	res.push_back(std::move(initialField));
	return res;
}
//
//static Sources buildResonantModeInitialField(
//	const FieldType& ft = E,
//	const Source::Polarization& p = Source::Polarization{ 0.0,0.0,1.0 },
//	const std::vector<std::size_t>& modes = { 1 },
//	const int dimension = 1)
//{
//	Sources res;
//	res.push_back(
//		std::move(std::make_unique<InitialField>(
//			SinusoidalMode{ dimension, modes }, ft, p) 
//		)
//	);
//	return res;
//}
//
//static Sources buildPlaneWave(
//	const Source::Polarization& p = Source::Polarization{ 0.0,0.0,1.0 },
//	const double spread = 0.1,
//	const double delay = 0.0,
//	const double normalization = 1.0,
//	const int dimension = 1)
//{
//	Sources res;
//	res.push_back(
//		std::move(std::make_unique<PlaneWave>(
//			TimeGaussian{ dimension, spread, normalization, delay }, p)
//		)
//	);
//	return res;
//}
//
//static Sources buildRightTravelingWaveInitialField(const Gaussian& gauss)
//{
//	Sources res;
//	res.push_back(std::move(std::make_unique<InitialField>(gauss, E, Source::Polarization{0.0, 1.0, 0.0})));
//	res.push_back(std::move(std::make_unique<InitialField>(gauss, H, Source::Polarization{0.0, 0.0, 1.0})));
//	return res;
//}

}
}
}

