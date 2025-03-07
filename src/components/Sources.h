#pragma once

#include "Types.h"
#include "math/Function.h"
#include "Spherical.h"
#include <functional>
#include "math/PhysicalConstants.h"
#include "math/Calculus.h"
#include <memory>
#include <cassert>
#include <float.h>

namespace maxwell {

class Source {
public:
	using Position = mfem::Vector;
	using Time = double;
	using Polarization = mfem::Vector;
	using Propagation = mfem::Vector;
	using CartesianAngles = mfem::Vector;

	virtual ~Source() = default;
	virtual std::unique_ptr<Source> clone() const = 0;

};

class InitialField : public Source {
public:
	InitialField(
		const Function&,
		const FieldType&,
		const Polarization&,
		const Position& center,
		const CartesianAngles& angles = CartesianAngles({ 0.0,0.0,0.0 })
	);
	InitialField(const InitialField&);

	std::unique_ptr<Source> clone() const;

	double eval(
		const Position&, const Time&, 
		const FieldType&, const Direction&) const;

	FieldType& fieldType() { return fieldType_; }
	Polarization& polarization() { return polarization_; }
	Function* magnitude() { return magnitude_.get(); }
	Position& center() { return center_; }

private:
	std::unique_ptr<Function> magnitude_;
	FieldType fieldType_{ E };
	Polarization polarization_;
	Position center_;
};

class TotalField : public Source {
public:
	TotalField(const EHFieldFunction&, const Polarization&, const Propagation&, const FieldType&);
	TotalField(const TotalField&);

	std::unique_ptr<Source> clone() const;

	VectorTF eval(
		const Position&, const Time&) const;

private: 
	std::unique_ptr<EHFieldFunction> function_;
	FieldType fieldtype_;
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
		{
			for (auto& v : rhs) {
				v_.push_back(std::move(v));
			}
		}
		return *this;
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

	Source* add(std::unique_ptr<Source>&& newV)
	{
		v_.push_back(std::move(newV));
		return v_.back().get();
	}

	Source* add(const Source& newV)
	{
		v_.push_back(newV.clone());
		return v_.back().get();
	}

	std::size_t size() const
	{
		return v_.size();
	}
private:
	std::vector<std::unique_ptr<Source>> v_;
};

}