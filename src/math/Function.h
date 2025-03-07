#pragma once

#include <mfem.hpp>
#include <math.h>
#include <cmath>
#include "math/PhysicalConstants.h"
#include "evolution/Fields.h"
#include "components/Spherical.h"
#include "math/Calculus.h"
#include "components/Types.h"

namespace maxwell {

using VectorTF = std::array<std::array<double, 3>, 2>;

using Position = mfem::Vector;
using Time = double;
using Polarization = mfem::Vector;
using Propagation = mfem::Vector;
using CartesianAngles = mfem::Vector;

class Function {
public:

	virtual ~Function() = default;

	virtual std::unique_ptr<Function> clone() const = 0;

	virtual double eval(const Position&) const = 0;
	
	virtual int dimension() const = 0;

};

class EHFieldFunction {
public:

	virtual ~EHFieldFunction() = default;

	virtual std::unique_ptr<EHFieldFunction> clone() const = 0;

	virtual VectorTF eval(const Position&, const Time&) const = 0;

};

class Gaussian : public Function {
public:
	/** 
	* A Gaussian function from R^{dim} to R. 
	* @param spread controls the width.
	* @param center_ must be of the same Size as dim.
	* @param dim between 1 and 3 for the dimension.
	*/
	Gaussian(
		double spread, 
		const Position center = Position({0.0}),
		int dim = 1
	) :
		spread_{ spread },
		center_{ center },
		dimension_{ dim }
	{
		assert(center_.Size() == dimension_);
	}

	int dimension() const { return dimension_; }

	std::unique_ptr<Function> clone() const {
		return std::make_unique<Gaussian>(*this);
	}

	double eval(const Position& pos) const {
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return 	exp(-pow(pos[X] - center_[X], 2) /
					(2.0 * pow(spread_, 2))
				);
		case 2:
			return 	exp(
					-(pow(pos[X] - center_[X], 2.0)
					+ pow(pos[Y] - center_[Y], 2.0)) /
					(2.0 * pow(spread_, 2.0))
				);
		case 3:
			return exp(
					-(pow(pos[X] - center_[X], 2.0)
					+ pow(pos[Y] - center_[Y], 2.0)
					+ pow(pos[Z] - center_[Z], 2.0)) /
					(2.0 * pow(spread_, 2.0))
				);
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:
	double spread_{ 2.0 };
	Position center_;
	int dimension_;
};

class SinusoidalMode : public Function {
public:
	SinusoidalMode(std::vector<std::size_t> modes) :
		modes_{ modes }
	{}

	int dimension() const { return (int) modes_.size(); }

	std::unique_ptr<Function> clone() const {
		return std::make_unique<SinusoidalMode>(*this);
	}

	double eval(const Position& pos) const
	{
		double res{ 1.0 };
		for (auto d{ 0 }; d < modes_.size(); ++d) {
			res *= sin(modes_[d] * M_PI * pos[d]);
		}
		return res;
	}

private:
	std::vector<std::size_t> modes_;
};

class BesselJ6 : public Function {
public:
	BesselJ6() {};

	std::unique_ptr<Function> clone() const {
		return std::make_unique<BesselJ6>(*this);
	}

	int dimension() const { return 2; }

	double eval(const Position& pos) const
	{
		double alpha6 = 13.589290170541217;
		double besselj6;
		#ifdef _WIN32 
			besselj6 = _jn(6, alpha6 * std::sqrt(std::pow(pos[0], 2.0) + std::pow(pos[1], 2.0)));
		#elif __linux__
			besselj6 = jn(6, alpha6 * std::sqrt(std::pow(pos[0], 2.0) + std::pow(pos[1], 2.0)));
		#endif
		double cosinetheta = std::cos(6.0 * std::atan2(pos[1], pos[0]));
		double cosinetime = 1.0; // t = 0 in cos(alpha6 * t)
		return besselj6 * cosinetheta * cosinetime;
	}

};

class DipoleDerivGauss : public EHFieldFunction {
public:
	DipoleDerivGauss(const double length, const double gaussianSpread, const double gaussdelay) :
		len_(length),
		spread_(gaussianSpread),
		gaussdelay_(gaussdelay)
	{
	}
	DipoleDerivGauss(const DipoleDerivGauss&);

	std::unique_ptr<EHFieldFunction> clone() const {
		return std::make_unique<DipoleDerivGauss>(*this);
	}

	VectorTF eval(
		const Position& p, const Time& t) const {

		VectorTF res;
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

		auto spreadterm = spread_ * std::sqrt(2.0);
		auto spreadsqrt2 = spreadterm * spreadterm;

		auto expArg = (t - gaussdelay_) / (spreadsqrt2);
		auto expArg2 = expArg * expArg;

		auto iret = -expArg * std::exp(-expArg / 2.0);
		auto diret = -iret * 2.0 * expArg / spreadsqrt2;
		auto doublediret = diret * (-2.0) * expArg / spreadsqrt2 + iret * (-2.0) / spreadsqrt2 / spreadsqrt2;

		const auto& ifpe0 = physicalConstants::invFourPiEps0_SI;
		const auto& ifp = physicalConstants::invFourPi;
		
		auto er = len_ * ifpe0 * 2.0 * cost * (iret / radius3 + diret / (cs * radius2));
		auto et = len_ * ifpe0 * sint * (iret / radius3 + diret / (cs * radius2) + doublediret / (cs2 * radius));
		auto hp = len_ * ifp * sint * (doublediret / (radius * cs) + diret / radius2);

		std::vector<double> resFieldE, resFieldH;

		resFieldE = pos.convertSphericalVectorFieldToCartesian(er, et, 0.0);
		res[E][0] = resFieldE[0];
		res[E][1] = resFieldE[1];
		res[E][2] = resFieldE[2];

		resFieldH = pos.convertSphericalVectorFieldToCartesian(0.0, 0.0, hp);
		res[H][0] = resFieldH[0];
		res[H][1] = resFieldH[1];
		res[H][2] = resFieldH[2];

		return res;

	}

private:
	double len_;
	double spread_;
	double gaussdelay_;

};

class Planewave : public EHFieldFunction {
public:
	Planewave(const Function& mag, const Polarization& p, const Propagation& dir, const FieldType& ft) :
		magnitude_{ mag.clone() },
		polarization_{ p },
		propagation_{ dir },
		fieldtype_{ ft }
	{}

	Planewave(const Planewave& rhs) :
		magnitude_{ rhs.magnitude_->clone() },
		polarization_{ rhs.polarization_ },
		propagation_{ rhs.propagation_ },
		fieldtype_{ rhs.fieldtype_ }
	{
		assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
		assert(std::abs(1.0 - propagation_.Norml2()) <= TOLERANCE);
	}

	std::unique_ptr<EHFieldFunction> clone() const {
		return std::make_unique<Planewave>(*this);
	};

	VectorTF eval(
		const Position& p, const Time& t) const {
		assert(p.Size() <= 3);

		VectorTF res;
		Polarization ePol(3), hPol(3);

		for (auto ft : { E, H }) {
			switch (fieldtype_) {
			case E:
				if (ft == E) {
					ePol = polarization_;
				}
				else {
					ePol = crossProduct(propagation_, polarization_);
				}
				break;
			case H:
				if (ft == H) {
					hPol = polarization_;
				}
				else {
					hPol = crossProduct(polarization_, propagation_);
				}
				break;
			}
		}

		Position pos(3);
		if (p.Size() == 1) {
			pos = Position({ p[0],  0.0, 0.0 });
		}
		else if (p.Size() == 2) {
			pos = Position({ p[0], p[1], 0.0 });
		}
		else {
			pos = p;
		}

		auto phaseDelay{ pos * propagation_ / physicalConstants::speedOfLight };

		for (auto d : { X, Y, Z }) {
			res[E][d] = magnitude_->eval(Position({ phaseDelay - t })) * ePol[d];
			res[H][d] = magnitude_->eval(Position({ phaseDelay - t })) * hPol[d];
		}

		return res;
	};

private:
	std::unique_ptr<Function> magnitude_;
	Polarization polarization_;
	Propagation propagation_;
	FieldType fieldtype_;
};

}
