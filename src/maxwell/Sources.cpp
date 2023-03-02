#include "Sources.h"

namespace maxwell {

InitialField::InitialField(
	const MathFunction& f, 
	const FieldType& fT, 
	const Polarization& p,
	const Position& centerIn) :
	function_{ f.clone() },
	fieldType{ fT },
	polarization{ p },
	center{ centerIn }
{}

InitialField::InitialField(const InitialField& rhs) :
	function_{ rhs.function_->clone() },
	fieldType{ rhs.fieldType },
	polarization{ rhs.polarization }
{}

std::unique_ptr<Source> InitialField::clone() const
{
	return std::make_unique<InitialField>(*this);
}

double InitialField::eval(const Position& p, Time t) const
{
	return function_->eval(p, t);
}

PlaneWave::PlaneWave(const MathFunction& f, const Polarization& p):
	function_{ f.clone() },
	polarization{ p }
{}

PlaneWave::PlaneWave(const PlaneWave& rhs) :
	function_{ rhs.function_->clone() },
	polarization{ rhs.polarization }
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