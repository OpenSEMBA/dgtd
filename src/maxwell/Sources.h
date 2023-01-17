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

	virtual double eval(const Position&, const Time& time) const = 0;

};

class InitialField : public Source {
public:
	InitialField(const MathFunction& f, const FieldType& fT, const Direction& d) :
		function_{ f.clone() },
		fieldType_{ fT },
		direction_{ d }
	{}

	std::unique_ptr<Source> clone() const {
		return std::make_unique<InitialField>(*this);
	}

	double eval(const Position&, const Time&) const;

private:
	FieldType fieldType_{ E };
	Direction direction_{ X };
	
	std::unique_ptr<MathFunction> function_;

	const void checkInputArguments();
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