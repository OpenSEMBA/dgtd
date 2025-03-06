#include "Sources.h"

#include <functional>
#include "math/PhysicalConstants.h"
#include "math/Calculus.h"
#include <memory>
#include <cassert>
#include <float.h>

#ifdef __linux__
#ifndef DBL_EPSILON
#	define DBL_EPSILON 2.2204460492503131e-16
#endif
#endif

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

Dipole::Dipole(const double length, const Gaussian& gaussian) :
	len_(length)
{
	gaussian_ = std::make_unique<Gaussian>(gaussian.clone());
}

std::unique_ptr<Source> Dipole::clone() const
{
	return std::make_unique<Dipole>(*this);
}

double Dipole::eval(
	const Position& p, const Time& t,
	const FieldType& f, const Direction& d) const
{

	SphericalVector pos(p);

	const auto& cs = physicalConstants::speedOfLight_SI;
	auto cs2 = cs * cs;
	auto delay = t - pos.radius / cs;

	auto radius = pos.radius;
	auto radius2 = radius * radius;
	auto radius3 = radius2 * radius;

	auto invSpeed = 1.0 / cs;
	auto invSpeed2 = invSpeed * invSpeed;
	auto invSpeedRho = invSpeed / radius;
	auto invSpeed2Rho = invSpeed2 / radius;
	auto invSpeedRho2 = invSpeed / radius2;

	auto sint = std::sin(pos.theta);
	auto cost = std::cos(pos.theta);

	const auto& spread = gaussian_->getSpread();
	auto spreadterm = spread * std::sqrt(2.0);
	auto spreadsqrt2 = spreadterm * spreadterm;

	auto expArg =  (t - delay) / (spreadsqrt2);
	auto expArg2 = expArg * expArg;

	auto iret = -expArg * std::exp(-expArg/2.0);
	auto diret = -iret * 2.0 * expArg / spreadsqrt2;
	auto doublediret = diret * (-2.0) * expArg / spreadsqrt2 + iret * (-2.0) / spreadsqrt2 / spreadsqrt2;

	const auto& ifpe0 = physicalConstants::invFourPiEps0_SI;
	const auto& ifp = physicalConstants::invFourPi;

	std::vector<double> resField;
	switch (f) {
	case FieldType::E:
		{
			auto er = len_ * ifpe0 * 2.0 * cost * (iret / radius3 + diret / (cs * radius2));
			auto et = len_ * ifpe0 * sint * (iret / radius3 + diret / (cs * radius2) + doublediret / (cs2 * radius));
			resField = pos.convertSphericalVectorFieldToCartesian(er, et, 0.0);
			switch (d) {
			case X:
				return resField[0];
			case Y:
				return resField[1];
			case Z:
				return resField[2];
			}
		}
	case FieldType::H: {
		auto hp = len_ * ifp * sint * (doublediret / (radius * cs) + diret / radius2);
		resField = pos.convertSphericalVectorFieldToCartesian(0.0, 0.0, hp);
		switch (d) {
		case X:
			return resField[0];
		case Y:
			return resField[1];
		case Z:
			return resField[2];
		}
		}
	}

	

}

}
