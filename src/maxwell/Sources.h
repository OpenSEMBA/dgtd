#pragma once

#include <functional>
#include "Types.h"
#include "Model.h"

namespace maxwell {

class Source {
public:
	Source(Model& model, double spread, double delay, Direction& d, FieldType& ft);

	FunctionCoefficient getFunction(double time, Direction d, FieldType ft);
	std::tuple<FieldType, Direction, double>& getGaussianFunctionTuple(const Position& pos);

private:
	double function_;
	double spread_;
	double delay_;
	Direction d_;
	FieldType ft_;
	double time_;
	Vector center_;
	Vector normalizedPos_;
	Vector minBB_, maxBB_;

};
}