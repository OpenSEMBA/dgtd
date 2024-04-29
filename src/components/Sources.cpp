#include "Sources.h"

#include <functional>
#include "math/PhysicalConstants.h"
#include "math/Calculus.h"
#include <memory>
#include <cassert>
#include <float.h>

#define DBL_EPSILON 2.2204460492503131e-016

namespace maxwell {

constexpr double TOLERANCE = 10.0 * DBL_EPSILON;

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
	assert(p.Size() == center_.Size());
	
	if (f != fieldType_) {
		return 0.0;
	}

	Position pos(p.Size());
	for (int i{ 0 }; i < p.Size(); ++i) {
		pos[i] = p[i] - center_[i];
	}
	
	return magnitude_->eval(pos) * polarization_[d];
}

Planewave::Planewave(
	const Function& mag, 
	const Polarization& p, 
	const Propagation& dir,
	const FieldType& ft):
	magnitude_{ mag.clone() },
	polarization_{ p },
	propagation_{ dir },
	fieldtype_{ ft }
{}

Planewave::Planewave(const Planewave& rhs) :
	magnitude_{ rhs.magnitude_->clone() },
	polarization_{ rhs.polarization_ },
	propagation_{ rhs.propagation_ },
	fieldtype_{ rhs.fieldtype_ }
{
	assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
	assert(std::abs(1.0 - propagation_.Norml2()) <= TOLERANCE);
}

std::unique_ptr<Source> Planewave::clone() const
{
	return std::make_unique<Planewave>(*this);
}

double Planewave::eval(
	const Position& p, const Time& t,
	const FieldType& f, const Direction& d) const
{
	assert(p.Size() <= 3);
	assert(f == E || f == H);
	assert(d == X || d == Y || d == Z);
	
	mfem::Vector fieldPol(3);

	switch (fieldtype_) {
	case E:
		if (f == E) {
			fieldPol = polarization_;
		}
		else {
			fieldPol = crossProduct(propagation_, polarization_);
		}
		break;
	case H:
		if (f == H) {
			fieldPol = polarization_;
		}
		else {
			fieldPol = crossProduct(polarization_, propagation_);
		}
		break;
	}

	mfem::Vector pos(3);
	if (p.Size() == 1) {
		pos = mfem::Vector({ p[0],  0.0, 0.0 });
	}
	else if (p.Size() == 2) {
		pos = mfem::Vector({ p[0], p[1], 0.0 });
	}
	else {
		pos = p;
	}
	auto phaseDelay{ pos * propagation_ / physicalConstants::speedOfLight };

	return magnitude_->eval(mfem::Vector({phaseDelay - t})) * fieldPol[d];
}

}