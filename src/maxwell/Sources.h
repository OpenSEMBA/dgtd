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

	virtual ~Source() = default;
	virtual std::unique_ptr<Source> clone() const = 0;

	virtual double eval(const Position&, Time) const = 0;

};

class InitialField : public Source {
public:
	InitialField(const MathFunction& f, const FieldType& fT, const Direction& d);
	InitialField(const InitialField& rhs);

	std::unique_ptr<Source> clone() const;

	double eval(const Position&, Time) const;

	FieldType fieldType{ E };
	Direction direction{ X };

private:
	std::unique_ptr<MathFunction> function_;
};

class PlaneWave : public Source {
public:
	PlaneWave(const MathFunction& f, const Direction& d);
	PlaneWave(const PlaneWave& rhs);

	std::unique_ptr<Source> clone() const;

	double eval(const Position&, Time) const;

	Direction direction{ X };

private: 
	std::unique_ptr<MathFunction> function_;
};

using Sources = std::vector<std::unique_ptr<Source>>;

}