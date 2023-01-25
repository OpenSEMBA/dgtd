#include "Sources.h"

namespace maxwell {

InitialField::InitialField(const MathFunction& f, const FieldType& fT, const Direction& d) :
	function_{ f.clone() },
	fieldType{ fT },
	direction{ d }
{}

InitialField::InitialField(const InitialField& rhs) :
	function_{ rhs.function_->clone() },
	fieldType{ rhs.fieldType },
	direction{ rhs.direction }
{}

std::unique_ptr<Source> InitialField::clone() const
{
	return std::make_unique<InitialField>(*this);
}

double InitialField::eval(const Position& p, Time t) const
{
	return function_->eval(p, t);
}

PlaneWave::PlaneWave(const MathFunction& f, const Direction& d):
	function_{ f.clone() },
	direction{ d }
{}

PlaneWave::PlaneWave(const PlaneWave& rhs) :
	function_{ rhs.function_->clone() },
	direction{ rhs.direction }
{}

std::unique_ptr<Source> PlaneWave::clone() const
{
	return std::make_unique<PlaneWave>(*this);
}

double PlaneWave::eval(const Position& p, Time t) const
{
	return function_->eval(p, t);
}

}