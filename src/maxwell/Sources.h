#pragma once

#include <functional>
#include <mfem.hpp>
#include "Types.h"

namespace maxwell {

class Source {
public:
	using Position = mfem::Vector;
	
	FieldType fieldType{E};
	Direction direction{X};
	InitialFieldType initialFT;
	Position center;

	virtual ~Source() = default;
	virtual std::unique_ptr<Source> clone() const = 0;

	virtual double eval3D(const mfem::Vector&) const = 0;
	virtual double eval2D(const mfem::Vector&) const = 0;
	virtual double eval1D(const mfem::Vector&) const = 0;

};

class GaussianInitialField : public Source {
public:
	using Position = mfem::Vector;

	std::function<double(const GaussianInitialField::Position&)> f = 0; 

	GaussianInitialField(
		const FieldType& ft, 
		const Direction& d,
		const double spread, 
		const double normalization,
		const Position center
	);

	std::unique_ptr<Source> clone() const {
		return std::make_unique<GaussianInitialField>(*this);
	}
	
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

	std::unique_ptr<Source> clone() const {
		return std::make_unique<PlanarSinusoidalInitialField>(*this);
	}

	double eval3D(const mfem::Vector&) const;
	double eval2D(const mfem::Vector&) const;
	double eval1D(const mfem::Vector&) const;

private:

	std::vector<std::size_t> modes_{ {0,0,0} };
	double coefficient_{ 1.0 };

	const void assembleModesVector(std::vector<std::size_t> modes);
};

using Sources = std::vector<std::unique_ptr<Source>>;

}