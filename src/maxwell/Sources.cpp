#include "Sources.h"
#include <math.h>

namespace maxwell {

Source::Source(Model& model, double spread, double delay, Direction& d, FieldType& ft) :
	spread_(spread),
	delay_(delay)
{
	model.getMesh().GetBoundingBox(minBB_, maxBB_, 0);
	for (int i = 0; i < minBB_.Size(); i++) {
		center_[i] = (maxBB_[i] - minBB_[i]) / 0.5;
	}
	normalizedPos_.SetSize(model.getMesh().Dimension());
	normalizedPos_ = 0.0;

	funcPack_.ft_ = ft;
	funcPack_.d_ = d;	
	funcPack_.fc_ = FunctionCoefficient(buildGaussianFunction);
	
}

double Source::buildGaussianFunction(const Position& pos)
{
	
	for (int i = 0; i < normalizedPos_.Size(); i++) {
		normalizedPos_[i] = 2 * (pos[i] - center_[i]) / (maxBB_[i] - minBB_[i]);
	}
	return (1.0 / (pow(spread_, 3.0) * pow(2.0 * M_PI, 3.0 / 2.0))) *
		exp(-1.0 * (pow(normalizedPos_[X], 2.0) + pow(normalizedPos_[Y], 2.0) + pow(normalizedPos_[Z], 2.0)) /
			(2.0 * pow(spread_, 2.0)));
}

Source::FunctionPackage& Source::buildGaussianFunctionPackage(FieldType& ft, Direction& d, std::function<double(const Position&)> f)
{
	FunctionCoefficient aux(f);
	FunctionPackage res(ft,d,aux);

	return res;
}


}