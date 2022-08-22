#pragma once

#include <functional>
#include <mfem.hpp>

#include "Types.h"

namespace maxwell {

class GaussianInitialField {
public:
	using Position = mfem::Vector;

	GaussianInitialField(
		const FieldType& ft, 
		const Direction& d, 
		const double spread, 
		const double normalization,
		const Position center
	);

	double evalGaussianFunction3D(const Position&) const;
	double evalGaussianFunction2D(const Position&) const;
	double evalGaussianFunction1D(const Position&) const;
	
	FieldType getFieldType() const { return fieldType_; }
	Direction getDirection() const { return direction_; }

private:
	FieldType fieldType_{E};
	Direction direction_{X};
	double spread_{2.0};
	double normalization_{1.0};
	Position center_;

	const void checkInputArguments();
};

using Sources = std::vector<GaussianInitialField>;

}