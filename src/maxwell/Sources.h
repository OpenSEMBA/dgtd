#pragma once

#include <functional>
#include <mfem.hpp>
#include "Types.h"

namespace maxwell {

class Source {
public:
	using Position = mfem::Vector;


	FieldType getFieldType() const { return fieldType_; }
	Direction getDirection() const { return direction_; }
	Position getCenter() const { return center_; }
	InitialFieldType getInitialFieldType() const { return initialFT_; }	

	double eval3D(const mfem::Vector&) const;
	double eval2D(const mfem::Vector&) const;
	double eval1D(const mfem::Vector&) const;

	const void setFieldType(const FieldType ft);
	const void setDirection(const Direction d);
	const void setCenter(const Position center);
	const void setInitialFieldType(const InitialFieldType ift);

private:

	FieldType fieldType_{E};
	Direction direction_{X};
	InitialFieldType initialFT_;
	Position center_;
};

class GaussianInitialField : public Source {
public:
	using Position = mfem::Vector;

	GaussianInitialField(
		const FieldType& ft, 
		const Direction& d,
		const double spread, 
		const double normalization,
		const Position center
	);
	
	double eval3D(const mfem::Vector&) const;
	double eval2D(const mfem::Vector&) const;
	double eval1D(const mfem::Vector&) const;

private:
	double spread_{2.0};
	double normalization_{1.0};

	const void checkInputArguments();
};

class PlanarSinusoidalInitialField : public Source {
public:
	using Position = mfem::Vector;

	PlanarSinusoidalInitialField(
		const FieldType& ft,
		const Direction& d,
		const std::vector<std::size_t> modes,
		const double coefficient,
		const Position center
	);

	double eval3D(const mfem::Vector&) const;
	double eval2D(const mfem::Vector&) const;
	double eval1D(const mfem::Vector&) const;

private:

	std::vector<std::size_t> modes_{ {0,0,0} };
	double coefficient_{ 1.0 };

	const void assembleModesVector(std::vector<std::size_t> modes);
};

using Sources = std::vector<Source>;

}