#pragma once

#include <mfem.hpp>

#include "Types.h"

#include <math.h>

namespace maxwell {

class MathFunction {
public:
	virtual ~MathFunction() = default;

	virtual std::unique_ptr<MathFunction> clone() const = 0;

	virtual double eval(const mfem::Vector&, double time = 0.0) const = 0;
	
	virtual int dimension() const = 0;

};

class Gaussian : public MathFunction {
public:
	Gaussian(double spatialSpread, int dim) :
		spatialSpread_{ spatialSpread },
		dimension_{ dim }
	{}

	int dimension() const 
	{ 
		return dimension_; 
	}

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<Gaussian>(*this);
	}

	double eval(const mfem::Vector& pos, double time = 0.0) const 
	{
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return 	exp(
					-pow(pos[X], 2) /
					(2.0 * pow(spatialSpread_, 2))
				);
		case 2:
			return 	exp(
					-(pow(pos[X], 2.0) 
					+ pow(pos[Y], 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 3:
			return exp(
					-(pow(pos[X], 2.0)
					+ pow(pos[Y], 2.0)
					+ pow(pos[Z], 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:
	double spatialSpread_{ 2.0 };
	int dimension_;
};

class TimeGaussian : public MathFunction {
public:
	TimeGaussian(double spatialSpread, double delay, int dimension) :
		spatialSpread_{ spatialSpread },
		dimension_{ dimension },
		delay_{ delay }
	{}

	int dimension() const { return dimension_; }

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<TimeGaussian>(*this);
	}

	double eval(const mfem::Vector& pos, double time = 0.0) const
	{
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return exp(
					-pow(pos[X] - (time - delay_), 2) /
					(2.0 * pow(spatialSpread_, 2))
				);
		case 2:
			return exp(
					-(pow(pos[X] - (time - delay_), 2.0)
					+ pow(pos[Y] - (time - delay_), 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 3:
			return exp(
					-(pow(pos[X] - (time - delay_), 2.0)
						+ pow(pos[Y] - (time - delay_), 2.0)
						+ pow(pos[Z] - (time - delay_), 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:
	int dimension_{ -1 };
	double spatialSpread_{ 2.0 };
	double delay_{ 0.0 };
};

class SinusoidalMode : public MathFunction {
public:
	SinusoidalMode(int dimension, std::vector<std::size_t> modes) :
		dimension_{ dimension },
		modes_{ modes }
	{}

	int dimension() const { return dimension_; }

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<SinusoidalMode>(*this);
	}

	double eval(const mfem::Vector& pos, double time) const 
	{
		double res{ 1.0 };
		for (auto d{ 0 }; d < dimension_; ++d) {
			res *= sin(modes_[d] * 2.0*M_PI * pos[d]);
		}
		return res;
	}

private:
	int dimension_{ -1 };
	std::vector<std::size_t> modes_;
};

}