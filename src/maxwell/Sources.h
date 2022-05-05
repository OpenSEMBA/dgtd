#pragma once

#include <functional>
#include "Types.h"
#include "Model.h"


namespace maxwell {

class Source {
public:

	Source(Model& model, double spread, double delay, Direction& d, FieldType& ft);

private:

	struct FunctionPackage {
		FieldType ft_;
		Direction d_;
		FunctionCoefficient fc_;

		FunctionPackage();
		FunctionPackage(FieldType& ft, Direction& d, FunctionCoefficient& fc);
	};

	double spread_;
	double delay_;
	Vector center_;
	Vector normalizedPos_;
	Vector minBB_, maxBB_;
	FunctionPackage funcPack_;

	double buildGaussianFunction(const Position& pos);
	FunctionPackage& buildGaussianFunctionPackage(FieldType&, Direction&, std::function<double(const Position&)> f);
};

}