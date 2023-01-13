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
	
	double eval3D(const Position&) const;
	double eval2D(const Position&) const;
	double eval1D(const Position&) const;

private:
	double spread_{2.0};
	double normalization_{1.0};

	const void checkInputArguments();
};

class SinusoidalInitialField : public Source {
public:
	using Position = mfem::Vector;

	SinusoidalInitialField(
		const FieldType& ft,
		const Direction& d,
		const std::vector<std::size_t> modes,
		const std::vector<double> coefficient,
		const Position center
	);

	std::unique_ptr<Source> clone() const {
		return std::make_unique<SinusoidalInitialField>(*this);
	}

	double eval3D(const Position&) const;
	double eval2D(const Position&) const;
	double eval1D(const Position&) const;

private:

	std::vector<std::size_t> modes_{ {0,0,0} };
	std::vector<double> coefficient_{ {1.0,1.0,1.0} };

	const void assembleModesVector(std::vector<std::size_t> modes);
};

class PlaneWave : public Source {
public:
	using Position = mfem::Vector;

	enum ExcitationType {
		Gaussian
	};

	PlaneWave(
		const Direction& d
	);

	std::unique_ptr<Source> clone() const {
		return std::make_unique<PlaneWave>(*this);
	}

	double eval3D(const Position&, double) const;
	double eval2D(const Position&, double) const;
	double eval1D(const Position&, double) const;

private: 

	int excitationType_ = Gaussian;
	double time_ = 0.0;

	double eval3D(const Position&) const; //These shouldn't be used, if not here there's issues with the virtuals from base. TODO
	double eval2D(const Position&) const; //These shouldn't be used, if not here there's issues with the virtuals from base. TODO
	double eval1D(const Position&) const; //These shouldn't be used, if not here there's issues with the virtuals from base. TODO

};

using Sources = std::vector<std::unique_ptr<Source>>;

}