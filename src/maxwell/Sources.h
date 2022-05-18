#pragma once

#include <functional>
#include "Types.h"
#include "Model.h"

namespace maxwell {


class Source {
public:
	Source(Model& model, const double spread, const double coeff, 
		const double devFromCenter, const Direction& d, const FieldType& ft);

	double evalGaussianFunction(const Position& pos) const;
	double evalGaussianFunction1D(const Position& pos) const;
	FieldType getFieldType() const { return fieldType_; }
	Direction getDirection() const { return direction_; }

private:

	FieldType fieldType_;
	Direction direction_;
	double spread_;
	double coeff_;
	double devFromCenter_;
	Vector minBB_, maxBB_;

	static Vector vectorAverage(const Vector& min, const Vector& max);
};

struct Sources {
public:

	void addSourceToVector(const Source& source) { sourceVector_.push_back(source); }
	std::vector<Source> getSourcesVector() const { return sourceVector_; }

	Sources() = default;

private:

	std::vector<Source> sourceVector_;

};

}