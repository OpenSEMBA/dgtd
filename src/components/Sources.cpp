#include "Sources.h"

#ifdef __linux__
#ifndef DBL_EPSILON
#	define DBL_EPSILON 2.2204460492503131e-16
#endif
#endif

namespace maxwell {

InitialField::InitialField(
	const Function& f, 
	const FieldType& fT, 
	const Polarization& p,
	const Position& centerIn,
	const CartesianAngles& angles) :
	function_{ f.clone() },
	fieldType_{ fT },
	polarization_{ p },
	center_{ centerIn }
{
	assert(std::abs(1.0 - polarization_.Norml2()) <= TOLERANCE);
}

InitialField::InitialField(const InitialField& rhs) :
	function_{ rhs.function_->clone() },
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
	Position center(p.Size());
	center = 0.0;
	for (int v{0}; v < center.Size(); v++){
		center[v] = center_[v];
	}
	for (int i{ 0 }; i < p.Size(); ++i) {
		pos[i] = p[i] - center[i];
	}
	
	return function_->eval(pos) * polarization_[d];
}

DeltaGapSource::DeltaGapSource(
	double magnitude,
	double spread,
	double t0,
	const mfem::Vector& polarization,
	bool derivative) :
	magnitude_{ magnitude },
	spread_{ spread },
	t0_{ t0 },
	derivative_{ derivative },
	polarization_{ polarization }
{
	if (!(magnitude_ > 0.0) || !std::isfinite(magnitude_)) {
		throw std::runtime_error("delta_gap magnitude must be > 0.");
	}
	if (!(spread_ > 0.0) || !std::isfinite(spread_)) {
		throw std::runtime_error("delta_gap spread must be > 0.");
	}
	if (!std::isfinite(t0_)) {
		throw std::runtime_error("delta_gap t0 must be finite.");
	}
	if (polarization_.Size() != 3 || !(polarization_.Norml2() > 0.0)) {
		throw std::runtime_error("delta_gap polarization must be a nonzero 3-vector.");
	}
	polarization_ /= polarization_.Norml2();
}

DeltaGapSource::DeltaGapSource(const DeltaGapSource& rhs) :
	magnitude_{ rhs.magnitude_ },
	spread_{ rhs.spread_ },
	t0_{ rhs.t0_ },
	derivative_{ rhs.derivative_ },
	polarization_{ rhs.polarization_ }
{}

std::unique_ptr<Source> DeltaGapSource::clone() const
{
	return std::make_unique<DeltaGapSource>(*this);
}

double DeltaGapSource::waveform(Time t) const
{
	return gaussianTimeSignal(t, t0_, spread_, magnitude_, derivative_);
}

double DeltaGapSource::eval(
	const Position&, const Time& t,
	const FieldType& ft, const Direction& d) const
{
	if (ft != E) {
		return 0.0;
	}
	if (d == X) {
		return waveform(t) * polarization_[0];
	}
	if (d == Y) {
		return waveform(t) * polarization_[1];
	}
	if (d == Z) {
		return waveform(t) * polarization_[2];
	}
	return 0.0;
}

TotalField::TotalField(
	const EHFieldFunction& func):
	function_{ func.clone() }
{}

TotalField::TotalField(const TotalField& rhs) :
	function_{ rhs.function_->clone() }
{
}

std::unique_ptr<Source> TotalField::clone() const
{
	return std::make_unique<TotalField>(*this);
}

double TotalField::eval(
	const Position& p, const Time& t,
	const FieldType& ft, const Direction& d) const
{
	return function_->eval(p, t, ft, d);
}

}
