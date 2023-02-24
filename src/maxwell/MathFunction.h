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

	std::vector<Direction> getFieldDirsByPolDir(const Direction& pol) const
	{
		std::vector<Direction> res(2);
		switch (pol)
		{
		case X:
			res[0] = Y, res[1] = Z;
			break;
		case Y:
			res[0] = Z, res[1] = X;
			break;
		case Z:
			res[0] = X, res[1] = Y;
			break;
		}
		return res;
	}

};

class Gaussian : public MathFunction {
public:
	Gaussian(int dimension, double spatialSpread, double normalization, const mfem::Vector& center, const Direction& polarization) :
		dimension_{ dimension },
		spatialSpread_{ spatialSpread },
		normalization_{ normalization },
		rotCenter_{ center },
		polarization_{polarization}
	{
	}

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<Gaussian>(*this);
	}

	double eval(const mfem::Vector& pos, double time = 0.0) const 
	{
		assert(dimension_ <= pos.Size());
		auto fieldDir{ getFieldDirsByPolDir(polarization_) };
		switch (dimension_) {
		case 1:
			return normalization_* 
				exp(
					-pow(pos[fieldDir[0]] - rotCenter_[fieldDir[0]], 2) /
					(2.0 * pow(spatialSpread_, 2))
				);
		case 2:
			return normalization_ * 
				exp(
					-(pow(pos[fieldDir[0]] - rotCenter_[fieldDir[0]], 2.0)
					+ pow(pos[fieldDir[1]] - rotCenter_[fieldDir[1]], 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 3:
			return normalization_ * 
				exp(
					-(pow(pos[fieldDir[0]] - rotCenter_[fieldDir[0]], 2.0)
					+ pow(pos[fieldDir[1]] - rotCenter_[fieldDir[1]], 2.0)
					+ pow(pos[fieldDir[2]] - rotCenter_[fieldDir[2]], 2.0)) /
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
	mfem::Vector rotCenter_;
	Direction polarization_{ Z };

};

class RotatedGaussian : public MathFunction {
public:
	RotatedGaussian(int dimension, double spatialSpread, double normalization, double angle, const mfem::Vector& center, const Direction& polarization, const Direction& rotation) :
		dimension_{ dimension },
		spatialSpread_{ spatialSpread },
		normalization_{ normalization },
		angle_{ angle },
		center_{ center },
		polarization_ { polarization },
		rotation_ { rotation }
	{}

	std::unique_ptr<MathFunction> clone() const {
		return std::make_unique<RotatedGaussian>(*this);
	}

	double eval(const mfem::Vector& pos, double time = 0.0) const
	{
		assert(dimension_ <= pos.Size());
		mfem::Vector rotCoords{ getRotatedCoords(pos) };
		auto fieldDir{ getFieldDirsByPolDir(polarization_) };
		switch (dimension_) {
		case 1:
			return normalization_ *
				exp(
					-pow((rotCoords[fieldDir[0]]), 2.0) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 2:
			return normalization_ *
				exp(
					-(pow(rotCoords[fieldDir[0]], 2.0)
					+ pow(rotCoords[fieldDir[1]], 2.0)) /
					(2.0 * pow(spatialSpread_, 2.0))
				);
		case 3:
			throw std::runtime_error("Rotated Gaussian Dim(3) - To Be Implemented.");
		default:
			throw std::runtime_error("Invalid dimension.");
		}
	}

private:

	int dimension_{ -1 };
	double spatialSpread_{ 2.0 };
	double normalization_{ 1.0 };
	double angle_{ 0.0 };
	mfem::Vector center_;
	Direction polarization_{ Z };
	Direction rotation_{ Z };
	mfem::DenseMatrix rotMat_;

	mfem::DenseMatrix getRotationMatrix()
	{
		switch (rotation_) {
		case X:
			return mfem::DenseMatrix({
				{1.0,         0.0,          0.0},
				{0.0, cos(angle_), -sin(angle_)},
				{0.0, sin(angle_),  cos(angle_)}});
		case Y:
			return mfem::DenseMatrix({
				{ cos(angle_), 0.0, sin(angle_)},
				{         0.0, 1.0,         0.0}, 
				{-sin(angle_), 0.0, cos(angle_)}});
		case Z:
			return mfem::DenseMatrix({
				{cos(angle_), -sin(angle_), 0.0},
				{sin(angle_),  cos(angle_), 0.0},
				{        0.0,          0.0, 1.0}});
		default:
			throw std::exception("Rotation Direction is not X, Y or Z.");
		}
	}

	mfem::Vector getRotatedCoords(const mfem::Vector& pos) const
	{
		mfem::Vector posSubCent(pos.Size()), res(pos.Size());
		for (int i{ 0 }; i < res.Size(); ++i) {
			posSubCent[i] = pos[i] - center_[i];
		}
		rotMat_.Mult(posSubCent, res);
		return res;
	}
		

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
	SinusoidalMode(int dimension, std::vector<std::size_t> modes) :
		dimension_{ dimension },
		modes_{ modes }
	{}

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