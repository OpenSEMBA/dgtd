#include "Sources.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace maxwell {

GaussianInitialField::GaussianInitialField(
	const FieldType& ft,
	const Direction& d,
	const double spread,
	const double normalization,
	const Position cnt) :
	spread_(spread),
	normalization_(normalization)
{
	fieldType = ft,
		direction = d,
		center = cnt,
		initialFT = InitialFieldType::Gaussian;
	checkInputArguments();
};

const void GaussianInitialField::checkInputArguments()
{
	if (spread_ < 0.0) {
		throw std::exception("Invalid spread value.");
	}
	if (normalization_ < 0.0) {
		throw std::exception("Invalid coeff value.");
	}
}

double GaussianInitialField::eval3D(const Position& pos) const
{
	return normalization_ * exp(-(pow(pos[X] - center[X], 2.0)
		+ pow(pos[Y] - center[Y], 2.0)
		+ pow(pos[Z] - center[Z], 2.0)) / (2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval2D(const Position& pos) const
{
	return normalization_ * exp(-(pow(pos[X] - center[X], 2.0)
		+ pow(pos[Y] - center[Y], 2.0)) / (2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval1D(const Position& pos) const
{
	return normalization_
		* exp(-pow(pos[X] - center[X], 2) / (2.0 * pow(spread_, 2)));
}

SinusoidalInitialField::SinusoidalInitialField(
	const FieldType& ft,
	const Direction& d,
	const std::vector<std::size_t> modes,
	const std::vector<double> coefficient,
	const Position cnt) :
	coefficient_(coefficient)
{
	fieldType = ft,
	direction = d,
	center = cnt,
	initialFT = InitialFieldType::PlanarSinusoidal;
	assembleModesVector(modes);
}

const void SinusoidalInitialField::assembleModesVector(std::vector<std::size_t> modes)
{
	for (int i = 0; i < modes.size(); i++) {
		modes_[i] = modes[i];
	}
}

double SinusoidalInitialField::eval3D(const Position& pos) const
{
	switch (direction) {
	case X:
		return sin(coefficient_[Y] * modes_[Y] * M_PI * pos[Y]) * sin(coefficient_[Z] * modes_[Z] * M_PI * pos[Z]);
	case Y:
		return sin(coefficient_[X] * modes_[X] * M_PI * pos[X]) * sin(coefficient_[Z] * modes_[Z] * M_PI * pos[Z]);
	case Z:
		return sin(coefficient_[X] * modes_[X] * M_PI * pos[X]) * sin(coefficient_[Y] * modes_[Y] * M_PI * pos[Y]);
	}
}

double SinusoidalInitialField::eval2D(const Position& pos) const
{
	return sin(coefficient_[X] * modes_[X] * M_PI * pos[X]) * sin(coefficient_[Y] * modes_[Y] * M_PI * pos[Y]);
}

double SinusoidalInitialField::eval1D(const Position& pos) const
{
	return sin(coefficient_[X] * modes_[X] * M_PI * pos[X]);
}
}