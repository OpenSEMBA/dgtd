#include "Sources.h"

#ifdef __linux__
#ifndef DBL_EPSILON
#	define DBL_EPSILON 2.2204460492503131e-16
#endif
#endif

namespace maxwell {

InitialField::InitialField(
	const Function& f, 
	const FieldType& fT, 
	const Polarization& p,
	const Position& centerIn,
	const CartesianAngles& angles) :
	function_{ f.clone() },
	fieldType_{ fT },
	polarization_{ p },
	center_{ centerIn }
{
	assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
}

InitialField::InitialField(const InitialField& rhs) :
	function_{ rhs.function_->clone() },
	fieldType_{ rhs.fieldType_ },
	polarization_{ rhs.polarization_ },
	center_{ rhs.center_ }
{}

std::unique_ptr<Source> InitialField::clone() const
{
	return std::make_unique<InitialField>(*this);
}

double InitialField::eval(
	const Position& p, const Time& t,
	const FieldType& f, const Direction& d) const
{
	
	if (f != fieldType_) {
		return 0.0;
	}

	Position pos(p.Size());
	for (int i{ 0 }; i < p.Size(); ++i) {
		pos[i] = p[i] - center_[i];
	}
	
	return function_->eval(pos) * polarization_[d];
}

TotalField::TotalField(
	const EHFieldFunction& func):
	function_{ func.clone() }
{}

TotalField::TotalField(const TotalField& rhs) :
	function_{ rhs.function_->clone() }
{
}

std::unique_ptr<Source> TotalField::clone() const
{
	return std::make_unique<TotalField>(*this);
}

double TotalField::eval(
	const Position& p, const Time& t,
	const FieldType& ft, const Direction& d) const
{
	return function_->eval(p, t, ft, d);
}

}
