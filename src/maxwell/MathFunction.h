#pragma once

#include <mfem.hpp>

#include "Types.h"

#include <math.h>

namespace maxwell {

class MathFunction {
public:
	virtual ~MathFunction() = default;

	virtual std::unique_ptr<MathFunction> clone() const = 0;

	virtual double eval(const mfem::Vector&) const = 0;
	
	virtual int dimension() const = 0;

};

class Gaussian : public MathFunction {
public:
	/** 
	* A Gaussian function from R^{dim} to R. 
	* @param spread controls the width.
	* @param center_ must be of the same Size as dim.
	* @param dim between 1 and 3 for the dimension.
	*/
	Gaussian(
		double spread, 
		const mfem::Vector center_ = mfem::Vector({0.0}), 
		int dim = 1
	) :
		spread_{ spread },
		center_{ center_ },
		dimension_{ dim }
	{
		assert(center_.Size() == dimension_);
	}

	int dimension() const { return dimension_; }

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<Gaussian>(*this);
	}

	double eval(const mfem::Vector& pos) const 
	{
		assert(dimension_ <= pos.Size());
		switch (dimension_) {
		case 1:
			return 	exp(-pow(pos[X], 2) /
					(2.0 * pow(spread_, 2))
				);
		case 2:
			return 	exp(
					-(pow(pos[X], 2.0) 
					+ pow(pos[Y], 2.0)) /
					(2.0 * pow(spread_, 2.0))
				);
		case 3:
			return exp(
					-(pow(pos[X], 2.0)
					+ pow(pos[Y], 2.0)
					+ pow(pos[Z], 2.0)) /
					(2.0 * pow(spread_, 2.0))
				);
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:
	double spread_{ 2.0 };
	int dimension_;
	mfem::Vector center_;
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

	double eval(const mfem::Vector& pos) const 
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