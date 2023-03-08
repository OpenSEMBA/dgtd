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
	using CartesianAngles = std::vector<double>;

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
		const CartesianAngles rotAngle = CartesianAngles({ 0.0,0.0,0.0 })
	);
	InitialField(const InitialField&);

	std::unique_ptr<Source> clone() const;

	double eval(const Position&, Time) const;

	FieldType fieldType{ E };
	Polarization polarization;
	Position center;
	CartesianAngles angles;

private:
	std::unique_ptr<MathFunction> function_;
};

class TimeVaryingField : public Source {
public:
	TimeVaryingField(const MathFunction& f, const Polarization& p, const Position& centerIn, const CartesianAngles& anglesIn);
	TimeVaryingField(const TimeVaryingField&);

	std::unique_ptr<Source> clone() const;

	double eval(const Position&, Time) const;

	Polarization polarization;
	Position center;
	CartesianAngles angles;

private: 
	std::unique_ptr<MathFunction> function_;
};

using Sources = std::vector<std::unique_ptr<Source>>;

}