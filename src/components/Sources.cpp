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
	magnitude_{ f.clone() },
	fieldType_{ fT },
	polarization_{ p },
	center_{ centerIn }
{
	assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
}

InitialField::InitialField(const InitialField& rhs) :
	magnitude_{ rhs.magnitude_->clone() },
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
	
	return magnitude_->eval(pos) * polarization_[d];
}

TotalField::TotalField(
	const EHFieldFunction& mag, 
	const Polarization& p, 
	const Propagation& dir,
	const FieldType& ft):
	function_{ mag.clone() },
	polarization_{ p },
	propagation_{ dir },
	fieldtype_{ ft }
{}

TotalField::TotalField(const TotalField& rhs) :
	function_{ rhs.function_->clone() },
	polarization_{ rhs.polarization_ },
	propagation_{ rhs.propagation_ },
	fieldtype_{ rhs.fieldtype_ }
{
	assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
	assert(std::abs(1.0 - propagation_.Norml2()) <= TOLERANCE);
}

std::unique_ptr<Source> TotalField::clone() const
{
	return std::make_unique<TotalField>(*this);
}

VectorTF TotalField::eval(
	const Position& p, const Time& t) const
{
	return function_->eval(p, t);
}

}
