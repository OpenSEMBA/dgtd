#pragma once

#include <functional>
#include <mfem.hpp>
#include "Types.h"
#include "MathFunction.h"

namespace maxwell {

class Source {
public:
	using Position = mfem::Vector;
	using Time = double;
	using Polarization = mfem::Vector;

	virtual ~Source() = default;
	virtual std::unique_ptr<Source> clone() const = 0;

	virtual double eval(const Position&, Time) const = 0;

};

class InitialField : public Source {
public:
	InitialField(
		const MathFunction&, 
		const FieldType&, 
		const Polarization&, 
		const Position& center,
		const double rotAngle = 0.0
	);
	InitialField(const InitialField&);

	std::unique_ptr<Source> clone() const;

	double eval(const Position&, Time) const;

	FieldType fieldType{ E };
	Polarization polarization;
	Position center;
	double rotAngle;

private:
	std::unique_ptr<MathFunction> function_;
};

class PlaneWave : public Source {
public:
	PlaneWave(const MathFunction&, const Polarization&);
	PlaneWave(const PlaneWave&);

	std::unique_ptr<Source> clone() const;

	double eval(const Position&, Time) const;

	Polarization polarization;

private: 
	std::unique_ptr<MathFunction> function_;
};

using Sources = std::vector<std::unique_ptr<Source>>;

}