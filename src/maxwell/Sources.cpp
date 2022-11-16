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

}