#include "Sources.h"
#include <math.h>

namespace maxwell {

Source::Source(Model& model, double spread, double delay, Direction& d, FieldType& ft) :
	spread_(spread),
	delay_(delay),
	d_(d),
	ft_(ft)
{
	model.getMesh().GetBoundingBox(minBB_, maxBB_, 0);
	for (int i = 0; i < minBB_.Size(); i++) {
		center_[i] = (maxBB_[i] - minBB_[i]) / 0.5;
	}
	normalizedPos_.SetSize(model.getMesh().Dimension());
	normalizedPos_ = 0.0;
}

std::tuple<FieldType, Direction, double>& Source::getGaussianFunctionTuple(const Position& pos)
{
	for (int i = 0; i < normalizedPos_.Size(); i++) {
		normalizedPos_[i] = 2 * (pos[i] - center_[i]) / (maxBB_[i] - minBB_[i]);
	}

	return std::make_tuple(ft_,d_,
		(1.0 / (pow(spread_, 3) * pow(2 * M_PI, 3.0 / 2.0))) * exp(-1.0 * (
		pow(normalizedPos_[X], 2) + pow(normalizedPos_[Y], 2) + pow(normalizedPos_[Z], 2)) /
		(2 * pow(spread_, 2))));
}

}