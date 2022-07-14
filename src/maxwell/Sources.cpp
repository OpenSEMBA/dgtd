#include "Sources.h"
#include <math.h>


namespace maxwell {



Source::Source(
	Model& model,
	const FieldType& ft,
	const Direction& d, 
	const double spread, 
	const double coeff, 
	const Vector devFromCenter) : 

	spread_(spread),
	coeff_(coeff),
	devFromCenter_(devFromCenter),
	fieldType_(ft),
	direction_(d)
{
	checkInputArguments(model);
};

const void Source::checkInputArguments(Model& model)
{
	if (spread_ < 0.0) {
		throw std::exception("Invalid spread value.");
	}
	if (coeff_ < 0.0) {
		throw std::exception("Invalid coeff value.");
	}
	if (model.getConstMesh().Dimension() != devFromCenter_.Size()) {
		throw std::exception("Mesh dimension and devFromCenter vector size are not the same.");
	}
	model.getMesh().GetBoundingBox(minBB_, maxBB_, 0);
	for (int i = 0; i < devFromCenter_.Size(); i++) {
		if (devFromCenter_[i] < minBB_[i] || devFromCenter_[i] > maxBB_[i]) {
			throw std::exception("Deviation from center cannot be smaller than min boundary or bigger than max boundary values.");
		}
	}
}

double Source::evalGaussianFunction3D(const Position& pos) const
{
	Vector center = vectorAverage(minBB_, maxBB_);
	Vector normalizedPos(pos.Size());
	normalizedPos = 0.0;
	for (int i = 0; i < normalizedPos.Size(); i++) {
		normalizedPos[i] = 2 * (pos[i] - center[i] - devFromCenter_[i]) / (maxBB_[i] - minBB_[i]);
	}
	return coeff_ * (1.0 / (pow(spread_, 2.0) * pow(2.0 * M_PI, 2.0 / 2.0))) *
		exp(-40 * (pow(normalizedPos[X], 2.0) + pow(normalizedPos[Y], 2.0) + pow(normalizedPos[Z], 2.0)) /
			(2.0 * pow(spread_, 2.0)));
}


double Source::evalGaussianFunction2D(const Position& pos) const
{
	Vector center = vectorAverage(minBB_, maxBB_);
	Vector normalizedPos(pos.Size());
	normalizedPos = 0.0;
	for (int i = 0; i < normalizedPos.Size(); i++) {
		normalizedPos[i] = 2 * (pos[i] - center[i] - devFromCenter_[i]) / (maxBB_[i] - minBB_[i]);
	}
	return coeff_ * (1.0 / (pow(spread_, 2.0) * pow(2.0 * M_PI, 2.0 / 2.0))) *
		exp(-40 * (pow(normalizedPos[X], 2.0) + pow(normalizedPos[Y], 2.0)) /
			(2.0 * pow(spread_, 2.0)));
}

double Source::evalGaussianFunction1D(const Position& pos) const
{
	double center = (minBB_[0] + maxBB_[0]) * 0.5 - devFromCenter_[0];
	double normalizedPos = 2 * (pos[0] - center) / (maxBB_[0] - minBB_[0]);
	return coeff_ * (1.0 / spread_ * sqrt(2.0 * M_PI)) *
		exp(-40 * pow(normalizedPos , 2.0) / pow(spread_, 2.0));

}
Vector Source::vectorAverage(const Vector& a, const Vector& b)
{
	Vector res = a;
	res.Add(1.0, b);
	res /= 2.0;
	return res;
}

}