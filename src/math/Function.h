#pragma once

#include <mfem.hpp>
#include <math.h>
#include <cmath>

namespace maxwell {

class Function {
public:
	virtual ~Function() = default;

	virtual std::unique_ptr<Function> clone() const = 0;

	virtual double eval(const mfem::Vector&) const = 0;
	
	virtual int dimension() const = 0;

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
		const mfem::Vector center = mfem::Vector({0.0}), 
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

	const double getSpread() const{ return spread_; }

	double eval(const mfem::Vector& pos) const 
	{
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
	mfem::Vector center_;
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

	double eval(const mfem::Vector& pos) const 
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

	double eval(const mfem::Vector& pos) const
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

}