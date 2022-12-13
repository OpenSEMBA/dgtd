#include "Sources.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace maxwell {

const void Source::setFieldType(const FieldType ft)
{
	fieldType_ = ft;
}

const void Source::setDirection(const Direction d)
{
	direction_ = d;
}

const void Source::setCenter(const Position center)
{
	center_ = center;
}

const void Source::setInitialFieldType(const InitialFieldType ift)
{
	initialFT_ = ift;
}

GaussianInitialField::GaussianInitialField(
	const FieldType& ft,
	const Direction& d, 
	const double spread, 
	const double normalization, 
	const Position center) : 
	spread_(spread),
	normalization_(normalization)
{
	setFieldType(ft);
	setDirection(d);
	setCenter(center);
	setInitialFieldType(InitialFieldType::Gaussian);

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
	return normalization_ * exp(-(pow(pos[X] - getCenter()[X], 2.0)
								+ pow(pos[Y] - getCenter()[Y], 2.0)
								+ pow(pos[Z] - getCenter()[Z], 2.0)) / (2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval2D(const Position& pos) const
{
	return normalization_ * exp( - (pow(pos[X] - getCenter()[X], 2.0)
								  + pow(pos[Y] - getCenter()[Y], 2.0)) /	(2.0 * pow(spread_, 2.0)));
}
double GaussianInitialField::eval1D(const Position& pos) const
{
	return normalization_ 
		* exp( - pow(pos[X] - getCenter()[X], 2) / (2.0*pow(spread_, 2)) );
}

PlanarSinusoidalInitialField::PlanarSinusoidalInitialField(
	const FieldType& ft,
	const Direction& d,
	const std::vector<std::size_t> modes,
	const double coefficient,
	const Position center) :
	coefficient_(coefficient)
{
	setFieldType(ft);
	setDirection(d);
	setCenter(center);
	setInitialFieldType(InitialFieldType::PlanarSinusoidal);
	assembleModesVector(modes);
}

const void PlanarSinusoidalInitialField::assembleModesVector(std::vector<std::size_t> modes)
{
	for (int i = 0; i < modes.size(); i++) {
		modes_[i] = modes[i];
	}
}

//double PlanarSinusoidalInitialField::eval3D(const Position& pos) const
//{
//	switch (direction_) {
//	case X:
//		return coefficient_ * sin(modes_[Y] * M_PI * pos[Y]) * sin(modes_[Z] * M_PI * pos[Z]);
//	case Y:
//		return coefficient_ * sin(modes_[X] * M_PI * pos[X]) * sin(modes_[Z] * M_PI * pos[Z]);
//	case Z:
//		return coefficient_ * sin(modes_[Y] * M_PI * pos[Y]) * sin(modes_[Y] * M_PI * pos[Y]);
//	}
//}
//
//double PlanarSinusoidalInitialField::eval2D(const Position& pos) const
//{
//	switch (direction_) {
//	case X:
//		return coefficient_ * sin(modes_[Y] * M_PI * pos[Y]);
//	case Y:
//		return coefficient_ * sin(modes_[X] * M_PI * pos[X]);
//	}
//}

double PlanarSinusoidalInitialField::eval1D(const Position& pos) const
{
	return coefficient_ * sin(modes_[X] * M_PI * pos[X]);
}
	

}