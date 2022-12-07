#include "Sources.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace maxwell {

GaussianInitialField::GaussianInitialField(
	const FieldType& ft,
	const Direction& d, 
	const double spread, 
	const double normalization, 
	const Position center) : 
	fieldType_(ft),
	direction_(d),
	spread_(spread),
	normalization_(normalization),
	center_(center)
{
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
	return normalization_ * exp(-(pow(pos[X] - center_[X], 2.0)
								+ pow(pos[Y] - center_[Y], 2.0) 
								+ pow(pos[Z] - center_[Z], 2.0)) / (2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval2D(const Position& pos) const
{
	return normalization_ * exp( - (pow(pos[X] - center_[X], 2.0)
								  + pow(pos[Y] - center_[Y], 2.0)) /	(2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval1D(const Position& pos) const
{
	return normalization_ 
		* exp( - pow(pos[X] - center_[X], 2) / (2.0*pow(spread_, 2)) );
}

PlanarSinusoidalInitialField::PlanarSinusoidalInitialField(
	const FieldType& ft,
	const Direction& d,
	const std::vector<int> modes,
	const Position center) :
	fieldType_(ft),
	direction_(d),
	modes_(modes),
	center_(center)
{
	checkInputArguments();
}

const void PlanarSinusoidalInitialField::checkInputArguments()
{
	if (modes_.size() != 3){
		throw std::exception("modes vector must have three entries, representing xmode, ymode and zmode.");
	}
}

double PlanarSinusoidalInitialField::eval3D(const Position& pos) const
{
	switch (direction_) {
	case X:
		return sin(modes_[Y] * M_PI * pos[Y]) * sin(modes_[Z] * M_PI * pos[Z]);
	case Y:
		return sin(modes_[X] * M_PI * pos[X]) * sin(modes_[Z] * M_PI * pos[Z]);
	case Z:
		return sin(modes_[Y] * M_PI * pos[Y]) * sin(modes_[Y] * M_PI * pos[Y]);
	}
}

double PlanarSinusoidalInitialField::eval2D(const Position& pos) const
{
	switch (direction_) {
	case X:
		return sin(modes_[Y] * M_PI * pos[Y]);
	case Y:
		return sin(modes_[X] * M_PI * pos[X]);
	}
}

double PlanarSinusoidalInitialField::eval1D(const Position& pos) const
{
	return sin(modes_[X] * M_PI * pos[X]);
}
	

}