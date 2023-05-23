#pragma once

#include <functional>
#include <mfem.hpp>
#include "Types.h"
#include "MathFunction.h"

namespace maxwell {

class Source {
public:
	using Position = mfem::Vector;
	using Time = double;
	using Polarization = mfem::Vector;
	using Propagation = mfem::Vector;
	using CartesianAngles = std::vector<double>;

	virtual ~Source() = default;
	virtual std::unique_ptr<Source> clone() const = 0;

	virtual double eval(
		const Position&, const Time&,
		const FieldType&, const Direction&) const = 0;
};

class InitialField : public Source {
public:
	InitialField(
		const MathFunction&,
		const FieldType&,
		const Polarization&,
		const Position& center_,
		const CartesianAngles rotAngle = CartesianAngles({ 0.0,0.0,0.0 })
	);
	InitialField(const InitialField&);

	std::unique_ptr<Source> clone() const;

	double eval(
		const Position&, const Time&, 
		const FieldType&, const Direction&) const;

private:
	std::unique_ptr<MathFunction> magnitude_;
	FieldType fieldType_{ E };
	Polarization polarization_;
	Position center_;
	CartesianAngles angles_;
};

class Planewave : public Source {
public:
	Planewave(const MathFunction&, const Polarization&, const Propagation&);
	Planewave(const Planewave&);

	std::unique_ptr<Source> clone() const;

	double eval(
		const Position&, const Time&,
		const FieldType&, const Direction&) const;

private: 
	std::unique_ptr<MathFunction> magnitude_;
	Polarization polarization_;
	Propagation propagation_;
};

class Sources {
public:
	Sources() = default;
	Sources(const Sources& rhs) 
	{
		for (auto& v : rhs) {
			v_.push_back(v->clone());
		}
	}
	Sources(Sources&& rhs)
	{
		for (auto& v: rhs) {
			v_.push_back(std::move(v));
		}
	}
	Sources& operator=(const Sources& rhs)
	{
		for (auto& v : rhs) {
			v_.push_back(v->clone());
		}
		return *this;
	}
	Sources& operator=(Sources&& rhs)
	{
		Sources res;
		{
			for (auto& v : rhs) {
				v_.push_back(std::move(v));
			}
		}
		return res;
	}

	std::vector<std::unique_ptr<Source>>::const_iterator begin() const
	{
		return v_.cbegin();
	}

	std::vector<std::unique_ptr<Source>>::iterator begin()
	{
		return v_.begin();
	}

	std::vector<std::unique_ptr<Source>>::iterator end()
	{
		return v_.end();
	}

	std::vector<std::unique_ptr<Source>>::const_iterator end() const
	{
		return v_.cend();
	}

	void add(std::unique_ptr<Source>&& newV)
	{
		v_.push_back(std::move(newV));
	}
private:
	std::vector<std::unique_ptr<Source>> v_;
};

}