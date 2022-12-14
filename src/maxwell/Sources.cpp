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
	return normalization_ * exp( - (pow(pos[X] - center[X], 2.0)
								  + pow(pos[Y] - center[Y], 2.0)) /	(2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval1D(const Position& pos) const
{
	return normalization_ 
		* exp( - pow(pos[X] - center[X], 2) / (2.0*pow(spread_, 2)) );
}

void GaussianInitialField::binder1D(GaussianInitialField& source) const
{
	f = std::bind(&GaussianInitialField::eval1D, source, std::placeholders::_1);
}

PlanarSinusoidalInitialField::PlanarSinusoidalInitialField(
	const FieldType& ft,
	const Direction& d,
	const std::vector<std::size_t> modes,
	const double coefficient,
	const Position cnt) :
	coefficient_(coefficient)
{
	fieldType = ft,
	direction = d,
	center = cnt,
	initialFT = InitialFieldType::PlanarSinusoidal;
	assembleModesVector(modes);
}

const void PlanarSinusoidalInitialField::assembleModesVector(std::vector<std::size_t> modes)
{
	for (int i = 0; i < modes.size(); i++) {
		modes_[i] = modes[i];
	}
}

double PlanarSinusoidalInitialField::eval3D(const Position& pos) const
{
	switch (direction) {
	case X:
		return coefficient_ * sin(modes_[Y] * M_PI * pos[Y]) * sin(modes_[Z] * M_PI * pos[Z]);
	case Y:
		return coefficient_ * sin(modes_[X] * M_PI * pos[X]) * sin(modes_[Z] * M_PI * pos[Z]);
	case Z:
		return coefficient_ * sin(modes_[Y] * M_PI * pos[Y]) * sin(modes_[Y] * M_PI * pos[Y]);
	}
}

double PlanarSinusoidalInitialField::eval2D(const Position& pos) const
{
	switch (direction) {
	case X:
		return coefficient_ * sin(modes_[Y] * M_PI * pos[Y]);
	case Y:
		return coefficient_ * sin(modes_[X] * M_PI * pos[X]);
	}
}

double PlanarSinusoidalInitialField::eval1D(const Position& pos) const
{
	return coefficient_ * sin(modes_[X] * M_PI * pos[X]);
}
	

}