#pragma once

#include <functional>
#include "Types.h"
#include "Model.h"

namespace maxwell {


class Source {
public:
	Source(Model& model, double spread, double delay, Direction& d, FieldType& ft);

	double evalGaussianFunction(const Position& pos) const;
	double evalGaussianFunction1D(const Position& pos) const;
	FieldType getFieldType() const { return fieldType_; }
	Direction getDirection() const { return direction_; }

private:

	FieldType fieldType_;
	Direction direction_;
	double spread_;
	double delay_;
	Vector minBB_, maxBB_;

	static Vector vectorAverage(const Vector& min, const Vector& max);
};

struct Sources {
public:

	void addSourceToVector(const Source& source) { sourceVector_.push_back(source); }
	std::vector<Source> getSourceVector() { return sourceVector_; }

	Sources() = default;

private:

	std::vector<Source> sourceVector_;

};

}