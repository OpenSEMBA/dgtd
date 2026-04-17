#pragma once

#include "math/PhysicalConstants.h"
#include "evolution/Fields.h"
#include "components/Spherical.h"
#include "math/Calculus.h"

#include <gsl/gsl_sf_bessel.h>  
#include <gsl/gsl_sf_legendre.h>

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

	virtual double eval(const Position&, const Time&, const FieldType&, const Direction&) const = 0;

};

class Gaussian : public Function {
public:
	/** 
	* A Gaussian function from R^{dim} to R. 
	* @param spread controls the width.
	* @param mean_ must be of the same Size as dim.
	* @param dim between 1 and 3 for the dimension.
	*/
	Gaussian(
		double spread, 
		const Position mean,
		int dim = 1
	) :
		spread_{ spread },
		mean_{ mean },
		dimension_{ dim }
	{
		assert(mean_.Size() == dimension_);
	}

	int dimension() const { return dimension_; }

	std::unique_ptr<Function> clone() const {
		return std::make_unique<Gaussian>(*this);
	}

	double eval(const Position& pos) const {
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return 	exp(-pow(pos[X] - mean_[X], 2) /
					(2.0 * pow(spread_, 2))
				);
		case 2:
			return 	exp(
					-(pow(pos[X] - mean_[X], 2.0)
					+ pow(pos[Y] - mean_[Y], 2.0)) /
					(2.0 * pow(spread_, 2.0))
				);
		case 3:
			return exp(
					-(pow(pos[X] - mean_[X], 2.0)
					+ pow(pos[Y] - mean_[Y], 2.0)
					+ pow(pos[Z] - mean_[Z], 2.0)) /
					(2.0 * pow(spread_, 2.0))
				);
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:
	double spread_;
	Position mean_;
	int dimension_;
};

/**
* A modulated Gaussian: exp(-((x-mean)^2)/(2*spread^2)) * cos(2*pi*freq*(x-mean))
* Uses a wide envelope (large spread) with a carrier frequency to place spectral
* content around freq without needing fine mesh resolution for a narrow Gaussian.
* freq is in normalized units (f_SI / c_SI).
*/
class ModulatedGaussian : public Function {
public:
	ModulatedGaussian(
		double spread,
		const Position mean,
		double freq,
		int dim = 1
	) :
		spread_{ spread },
		mean_{ mean },
		freq_{ freq },
		dimension_{ dim }
	{
		assert(mean_.Size() == dimension_);
	}

	int dimension() const { return dimension_; }

	std::unique_ptr<Function> clone() const {
		return std::make_unique<ModulatedGaussian>(*this);
	}

	double eval(const Position& pos) const {
		assert(dimension_ <= pos.Size());
		double arg = pos[X] - mean_[X];
		double envelope = exp(-arg * arg / (2.0 * spread_ * spread_));
		double carrier = cos(2.0 * M_PI * freq_ * arg);
		return envelope * carrier;
	}

private:
	double spread_;
	Position mean_;
	double freq_;
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

class SphericalBesselJ6 : public Function {
public:
	SphericalBesselJ6() {}

	std::unique_ptr<Function> clone() const override {
		return std::make_unique<SphericalBesselJ6>(*this);
	}

	int dimension() const override { return 3; }

	double eval(const Position& pos) const override
	{
		const double alpha6 = 13.589290170541217;
		const double x = pos[0];
		const double y = pos[1];
		const double z = pos[2];

		double r = pos.Norml2();
		if (r == 0.0) return 0.0;

		double theta = std::acos(z / r);
		double phi = std::atan2(y, x);

		double j6 = gsl_sf_bessel_jl(6, alpha6 * r);

		double P66 = gsl_sf_legendre_sphPlm(6, 6, std::cos(theta)); 
		double Y66_real = P66 * std::cos(6 * phi);

		double cosAtT0 = 1.0;

		return j6 * Y66_real * cosAtT0;
	}
};

class DerivGaussDipole : public EHFieldFunction {
public:
	DerivGaussDipole(const double length, const double gaussianSpread, const double gaussMean) :
		len_(length),
		gaussSpread_(gaussianSpread),
		gaussMean_(gaussMean)
	{
	}

	DerivGaussDipole(const DerivGaussDipole& rhs) :
		len_{ rhs.len_ },
		gaussSpread_{ rhs.gaussSpread_ },
		gaussMean_{ rhs.gaussMean_ }
	{
	}

	std::unique_ptr<EHFieldFunction> clone() const {
		return std::make_unique<DerivGaussDipole>(*this);
	}

	double eval(
		const Position& p, const Time& t,
		const FieldType& ft, const Direction& d) const {

		SphericalVector pos(p);

		const auto& cs = physicalConstants::speedOfLight;
		auto cs2 = cs * cs;

		auto radius = pos.radius / cs;
		auto radius2 = radius * radius;
		auto radius3 = radius2 * radius;

		auto invSpeed = 1.0 / cs;
		auto invSpeed2 = invSpeed * invSpeed;
		auto invSpeedRho = invSpeed / radius;
		auto invSpeed2Rho = invSpeed2 / radius;
		auto invSpeedRho2 = invSpeed / radius2;

		auto sint = std::sin(pos.theta);
		auto cost = std::cos(pos.theta);

		auto spreadsqrt2 = gaussSpread_ * std::sqrt(2.0) ;

		auto expArg = (t - gaussMean_ - pos.radius / cs) / spreadsqrt2;
		auto expArg2 = expArg * expArg;

		auto maxMagnitude = 1.0;
		auto scalingFactor = (spreadsqrt2 * std::exp(1.0) / 2.0) * maxMagnitude;

		auto iret = scalingFactor * std::exp(-expArg2);
		auto diret = -iret * 2.0 * expArg / spreadsqrt2;
		auto doublediret = diret * (-2.0) * expArg / spreadsqrt2 + iret * (-2.0) / spreadsqrt2 / spreadsqrt2;

		const auto& ifpe0 = physicalConstants::invFourPiEps0;
		const auto& ifp = physicalConstants::invFourPi;

		//auto er = 0.0;
		//auto et = (len_ / cs) * ifpe0 * sint * (doublediret / (cs2 * radius));
		//auto hp = (len_ / cs) * ifp * sint * (doublediret / (radius * cs));

		auto er = (len_ / cs) * ifpe0 * 2.0 * cost * (iret / radius3 + diret / (cs * radius2));
		auto et = (len_ / cs) * ifpe0 * sint * (iret / radius3 + diret / (cs * radius2) + doublediret / (cs2 * radius));
		auto hp = (len_ / cs) * ifp * sint * (doublediret / (radius * cs) + diret / radius2);

		std::vector<double> resField;
		switch (ft) {
			case E:
				resField = pos.convertSphericalVectorFieldToCartesian(er, et, 0.0);
				switch (d) {
					case X:
						return resField[0];
					case Y:
						return resField[1];
					case Z:
						return resField[2];
				}
			case H:
				resField = pos.convertSphericalVectorFieldToCartesian(0.0, 0.0, hp);
				switch (d) {
					case X:
						return resField[0];
					case Y:
						return resField[1];
					case Z:
						return resField[2];
				}
			default:
				throw std::runtime_error("Unknown FieldType as argument in DerivGaussDipole.");
		}
	}

private:
	double len_;
	double gaussSpread_;
	double gaussMean_;

};

class Planewave : public EHFieldFunction {
public:
	Planewave(const Function& func, const Polarization& p, const Propagation& dir, const FieldType& ft) :
		function_{ func.clone() },
		fieldtype_{ ft }
	{
		polarization_.SetSize(p.Size());
		propagation_.SetSize(dir.Size());
		for (auto d{ 0 }; d < polarization_.Size(); d++) {
			polarization_[d] = p[d] / p.Norml2();
			propagation_[d] = dir[d] / dir.Norml2();
		}
	}

	Planewave(const Planewave& rhs) :
		function_{ rhs.function_->clone() },
		polarization_{ rhs.polarization_ },
		propagation_{ rhs.propagation_ },
		fieldtype_{ rhs.fieldtype_ }
	{
		if (!(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE)) {
			throw std::runtime_error("Polarization is not normalised.");
		}
		if (!(std::abs(1.0 - propagation_.Norml2()) <= TOLERANCE)) {
			throw std::runtime_error("Propagation is not normalised.");
		}
	}

	std::unique_ptr<EHFieldFunction> clone() const {
		return std::make_unique<Planewave>(*this);
	};

	double eval(
		const Position& p, const Time& t,
		const FieldType& ft, const Direction& d) const {
		assert(p.Size() <= 3);

		double polDir;

		switch (fieldtype_) {
		case E:
			if (ft == E) {
				polDir = polarization_[d];
			}
			else {
				{
					const auto cross = crossProduct(propagation_, polarization_);
					polDir = cross[d];
				}
			}
			break;
		case H:
			if (ft == H) {
				polDir = polarization_[d];
			}
			else {
				{
					const auto cross = crossProduct(polarization_, propagation_);
					polDir = cross[d];
				}
			}
			break;
		}

		Position pos;

		if (p.Size() == 1) {
			pos = Position({ p[0],  0.0, 0.0 });
		}
		else if (p.Size() == 2) {
			pos = Position({ p[0], p[1], 0.0 });
		}
		else {
			pos.SetDataAndSize(p.GetData(), 3);
		}

		auto phaseDelay{ pos * propagation_ / physicalConstants::speedOfLight };

		return function_->eval(Position({ phaseDelay - t })) * polDir;

	};

private:
	std::unique_ptr<Function> function_;
	Polarization polarization_;
	Propagation propagation_;
	FieldType fieldtype_;
};


class TimeFunction {
public:

	virtual ~TimeFunction() = default;

	virtual std::unique_ptr<TimeFunction> clone() const = 0;

	virtual double eval(const Position&, const Time&) const = 0;

};

class TM2DSinusoidalMode : public TimeFunction
{
public:
	TM2DSinusoidalMode(const std::vector<std::size_t>& modes, const std::vector<double>& box_size)
    {
        int dim = modes.size();
        if(modes_.size() != box_size_.size()){
            modes_.size() > box_size_.size() ? dim = modes_.size() : dim = box_size_.size();
        }
        modes_.resize(dim);
        box_size_.resize(dim);
        for (auto d = 0; d < dim; d++){
            modes_[d] = modes[d];
            box_size_[d] = box_size[d];
        }
    }

	std::unique_ptr<TimeFunction> clone() const {
		return std::make_unique<TM2DSinusoidalMode>(*this);
	}

	double eval(const Position& pos, const Time& t) const
	{
		double w = M_PI;
		double root_factor = 0.0;
        for (auto d = 0; d < modes_.size(); d++){
            root_factor += std::pow(modes_[d] / box_size_[d], 2);
        }
		w *= std::sqrt(root_factor);

        double res = std::cos(w * t);
        for (auto d = 0; d < modes_.size(); d++){
            res *= std::sin(modes_[d] * M_PI * pos[d] / box_size_[d]);
        }

        return res;
	}

private:
	std::vector<std::size_t> modes_;
    std::vector<double> box_size_;
};

class TM2D_Dx_SinusoidalMode : public TimeFunction
{
public:
	TM2D_Dx_SinusoidalMode(const std::vector<std::size_t>& modes, const std::vector<double>& box_size)
    {
        int dim = modes.size();
        if(modes_.size() != box_size_.size()){
            modes_.size() > box_size_.size() ? dim = modes_.size() : dim = box_size_.size();
        }
        modes_.resize(dim);
        box_size_.resize(dim);
        for (auto d = 0; d < dim; d++){
            modes_[d] = modes[d];
            box_size_[d] = box_size[d];
        }
    }

	std::unique_ptr<TimeFunction> clone() const {
		return std::make_unique<TM2D_Dx_SinusoidalMode>(*this);
	}

	double eval(const Position& pos, const Time& t) const
	{
		double w = M_PI;
		double root_factor = 0.0;
        for (auto d = 0; d < modes_.size(); d++){
            root_factor += std::pow(modes_[d] / box_size_[d], 2);
        }
		w *= std::sqrt(root_factor);

		double constant_factor = modes_[0] * M_PI / (w * box_size_[0]);
        double res = constant_factor * std::cos(w * t) 
		* std::cos(modes_[0] * M_PI * pos[0] / box_size_[0])
		* std::sin(modes_[1] * M_PI * pos[1] / box_size_[1]);

        return res;
	}

private:
	std::vector<std::size_t> modes_;
    std::vector<double> box_size_;
};

class TM2D_Dy_SinusoidalMode : public TimeFunction
{
public:
	TM2D_Dy_SinusoidalMode(const std::vector<std::size_t>& modes, const std::vector<double>& box_size)
    {
        int dim = modes.size();
        if(modes_.size() != box_size_.size()){
            modes_.size() > box_size_.size() ? dim = modes_.size() : dim = box_size_.size();
        }
        modes_.resize(dim);
        box_size_.resize(dim);
        for (auto d = 0; d < dim; d++){
            modes_[d] = modes[d];
            box_size_[d] = box_size[d];
        }
    }

	std::unique_ptr<TimeFunction> clone() const {
		return std::make_unique<TM2D_Dy_SinusoidalMode>(*this);
	}

	double eval(const Position& pos, const Time& t) const
	{
		double w = M_PI;
		double root_factor = 0.0;
        for (auto d = 0; d < modes_.size(); d++){
            root_factor += std::pow(modes_[d] / box_size_[d], 2);
        }
		w *= std::sqrt(root_factor);

		double constant_factor = modes_[1] * M_PI / (w * box_size_[1]);
        double res = constant_factor * std::cos(w * t) 
		* std::sin(modes_[0] * M_PI * pos[0] / box_size_[0])
		* std::cos(modes_[1] * M_PI * pos[1] / box_size_[1]);

        return res;
	}

private:
	std::vector<std::size_t> modes_;
    std::vector<double> box_size_;
};

}
