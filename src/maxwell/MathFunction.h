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
};

class Gaussian : public MathFunction {
public:
	Gaussian(int dimension, double spatialSpread, double normalization, const mfem::Vector& center) :
		dimension_{ dimension },
		spatialSpread_{ spatialSpread },
		normalization_{ normalization },
		center_{ center }
	{}

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<Gaussian>(*this);
	}

	double eval(const mfem::Vector& pos, double time = 0.0) const 
	{
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return normalization_* 
				exp(
					-pow(pos[X] - center_[X], 2) / 
					(2.0 * pow(spatialSpread_, 2))
				);
		case 2:
			return normalization_ * 
				exp(
					-(pow(pos[X] - center_[X], 2.0)
					+ pow(pos[Y] - center_[Y], 2.0)) / 
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 3:
			return normalization_ * 
				exp(
					-(pow(pos[X] - center_[X], 2.0)
					+ pow(pos[Y] - center_[Y], 2.0)
					+ pow(pos[Z] - center_[Z], 2.0)) / 
					(2.0 * pow(spatialSpread_, 2.0))
				);
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:
	int dimension_{ -1 };
	double spatialSpread_{ 2.0 };
	double normalization_{ 1.0 };
	mfem::Vector center_;

};

class TimeGaussian : public MathFunction {
public:
	TimeGaussian(int dimension, double spatialSpread, double normalization, double delay) :
		dimension_{ dimension },
		spatialSpread_{ spatialSpread },
		normalization_{ normalization },
		delay_{ delay }
	{}

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<TimeGaussian>(*this);
	}

	double eval(const mfem::Vector& pos, double time = 0.0) const
	{
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return normalization_ *
				exp(
					-pow(pos[X] - (time - delay_), 2) /
					(2.0 * pow(spatialSpread_, 2))
				);
		case 2:
			return normalization_ *
				exp(
					-(pow(pos[X] - (time - delay_), 2.0)
						+ pow(pos[Y] - (time - delay_), 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 3:
			return normalization_ *
				exp(
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
	double normalization_{ 1.0 };
	double delay_{ 0.0 };

};

class SinusoidalMode : public MathFunction {
public:

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<SinusoidalMode>(*this);
	}

	double eval(const mfem::Vector& pos, double time) const 
	{
		double res{ 0.0 };
		for (auto d{ 0 }; d < dimension_; ++d) {
			res *= sin(modes_[d] * M_PI * pos[d]);
		}
		return res;
	}

private:
	int dimension_{ -1 };
	std::vector<std::size_t> modes_;
};

}