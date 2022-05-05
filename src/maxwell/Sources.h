#pragma once

#include <functional>
#include "Types.h"
#include "Model.h"


namespace maxwell {

struct ModelParam {
	Vector center_;
	Vector normalizedPos_;
	Vector minBB_, maxBB_;

	ModelParam();
	ModelParam(Model&);
};

struct FunctionPackage {
	FieldType ft_;
	Direction d_;
	FunctionCoefficient fc_;

	FunctionPackage();
	FunctionPackage(FieldType& ft, Direction& d, FunctionCoefficient& fc);
};


class Source {
public:

	Source(Model& model, double spread, double delay, Direction& d, FieldType& ft);

	FunctionPackage& getFunctionPackage() { return funcPack_; }
	ModelParam& getModelParam() { return modelParam_; }

private:

	double spread_;
	double delay_;
	FunctionPackage funcPack_;
	ModelParam modelParam_;

	double buildGaussianFunction(const Position& pos);
	FunctionPackage& buildGaussianFunctionPackage(FieldType&, Direction&, std::function<double(const Position&)> f);
};

}