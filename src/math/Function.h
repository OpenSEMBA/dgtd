#pragma once

#include "math/PhysicalConstants.h"
#include "evolution/Fields.h"
#include "components/Spherical.h"
#include "math/Calculus.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

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

	/// Evaluate 1D profile at normalized region coordinate xi in [0, 1],
	/// stretching the profile's native domain across the target mesh region.
	virtual double evalAtUnitCoord(double xi) const
	{
		Position p(1);
		p[0] = xi;
		return eval(p);
	}
	
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

	double spread() const { return spread_; }
	const Position& mean() const { return mean_; }

private:
	double spread_;
	Position mean_;
	int dimension_;
};

class ConstantProfile : public Function {
public:
	explicit ConstantProfile(double value, int dim = 1)
		: value_(value), dimension_(dim)
	{
	}

	int dimension() const override { return dimension_; }

	std::unique_ptr<Function> clone() const override
	{
		return std::make_unique<ConstantProfile>(*this);
	}

	double eval(const Position&) const override { return value_; }

	double evalAtUnitCoord(double) const override { return value_; }

private:
	double value_;
	int dimension_;
};

class GaussianProfile : public Function {
public:
	static constexpr double kDefaultCutoffFactor = 1.0e-2;

	GaussianProfile(
		double peak,
		double spread,
		const Position& center,
		double offset = 0.0,
		double cutoff_factor = kDefaultCutoffFactor)
		: peak_(peak),
		  offset_(offset),
		  cutoff_factor_(cutoff_factor),
		  z_cut_(computeZCut(cutoff_factor)),
		  gaussian_(spread, center, center.Size())
	{
	}

	int dimension() const override { return gaussian_.dimension(); }

	std::unique_ptr<Function> clone() const override
	{
		return std::make_unique<GaussianProfile>(*this);
	}

	double eval(const Position& pos) const override
	{
		return std::max(0.0, offset_ + peak_ * gaussian_.eval(pos));
	}

	double evalAtUnitCoord(double xi) const override
	{
		const double z = (2.0 * xi - 1.0) * z_cut_;
		return std::max(0.0, offset_ + peak_ * std::exp(-0.5 * z * z));
	}

	double cutoffFactor() const { return cutoff_factor_; }
	double zCut() const { return z_cut_; }

private:
	static double computeZCut(double cutoff_factor)
	{
		if (!(cutoff_factor > 0.0 && cutoff_factor < 1.0)) {
			throw std::runtime_error("GaussianProfile cutoff_factor must be in (0, 1).");
		}
		return std::sqrt(-2.0 * std::log(cutoff_factor));
	}

	double peak_;
	double offset_;
	double cutoff_factor_;
	double z_cut_;
	Gaussian gaussian_;
};

class PositiveGaussianDerivativeProfile : public Function {
public:
	static constexpr double kDefaultCutoffFactor = 1.0e-2;

	PositiveGaussianDerivativeProfile(
		double peak,
		double spread,
		const Position& center,
		int axis = 0,
		double offset = 0.0,
		double cutoff_factor = kDefaultCutoffFactor)
		: peak_(peak),
		  spread_(spread),
		  center_(center),
		  axis_(axis),
		  offset_(offset),
		  cutoff_factor_(cutoff_factor),
		  z_cut_(computeDerivativeZCut(cutoff_factor)),
		  gaussian_(spread, center, center.Size())
	{
	}

	int dimension() const override { return gaussian_.dimension(); }

	std::unique_ptr<Function> clone() const override
	{
		return std::make_unique<PositiveGaussianDerivativeProfile>(*this);
	}

	double eval(const Position& pos) const override
	{
		assert(axis_ >= 0 && axis_ < pos.Size());
		const double z = (pos[axis_] - center_[axis_]) / spread_;
		const double env = std::abs(z) * std::exp(0.5 * (1.0 - z * z));
		return std::max(0.0, offset_ + peak_ * env);
	}

	double evalAtUnitCoord(double xi) const override
	{
		const double z = (2.0 * xi - 1.0) * z_cut_;
		const double env = std::abs(z) * std::exp(0.5 * (1.0 - z * z));
		return std::max(0.0, offset_ + peak_ * env);
	}

	double cutoffFactor() const { return cutoff_factor_; }
	double zCut() const { return z_cut_; }

private:
	/// Solves z * exp(0.5 * (1 - z^2)) = cutoff_factor for z > 1 so that both
	/// derivative lobes peak at 1 (at z = +/-1) and decay to cutoff_factor at z = +/-z_cut.
	static double computeDerivativeZCut(double cutoff_factor)
	{
		if (!(cutoff_factor > 0.0 && cutoff_factor < 1.0)) {
			throw std::runtime_error("PositiveGaussianDerivativeProfile cutoff_factor must be in (0, 1).");
		}
		const double target = std::log(cutoff_factor);
		double z = std::sqrt(1.0 - 2.0 * target);
		for (int iter = 0; iter < 10; ++iter) {
			const double g = std::log(z) + 0.5 * (1.0 - z * z) - target;
			const double dg = (1.0 / z) - z;
			z -= g / dg;
		}
		return z;
	}

	double peak_;
	double spread_;
	Position center_;
	int axis_;
	double offset_;
	double cutoff_factor_;
	double z_cut_;
	Gaussian gaussian_;
};

// Evaluates RAM-C II 61 km stagnation-line omega_p(x) [rad/s].
// Expects x in metres with the nose at x = 0 and shock around x = -0.0073 m.
inline double ramc61kmStagnationOmegaP(const mfem::Vector& p)
{
	const double x = p(0);
	if (x > 0.0) {
		return 0.0;
	}
	constexpr double wp_inf = 9.697700e+08;
	constexpr double x_sh   = -7.309150e-03;
	constexpr double d_sh   =  3.211841e-04;
	constexpr double a0     =  6.865065e+11;
	constexpr double a1     = -4.836413e+13;
	constexpr double a2     = -2.041384e+15;
	constexpr double a3     =  8.678069e+17;
	constexpr double d_bl   =  4.296268e-04;
	constexpr double A_bl   =  2.639789e+11;

	const double step = 0.5 * (1.0 + std::tanh((x - x_sh) / d_sh));
	const double s    = x + 0.0066;
	const double env  = a0 + s * (a1 + s * (a2 + s * a3)) - A_bl * std::exp(x / d_bl);
	return wp_inf + step * std::max(0.0, env);
}

// Evaluates RAM-C II 61 km stagnation-line nu_c(x) [1/s] (Gupta closure).
inline double ramc61kmStagnationNuC_Gupta(const mfem::Vector& p)
{
	const double x = p(0);
	if (x > 0.0) {
		return 0.0;
	}
	constexpr double nu_inf =  2.732100e+07;
	constexpr double x_sh   = -7.482650e-03;
	constexpr double d_sh   =  4.340546e-04;
	constexpr double a0     =  5.765874e+09;
	constexpr double a1     =  6.977597e+11;
	constexpr double a2     = -1.255220e+14;
	constexpr double a3     =  8.677992e+15;
	constexpr double d_bl   =  2.616807e-04;
	constexpr double A_bl   =  2.012032e+09;

	const double step = 0.5 * (1.0 + std::tanh((x - x_sh) / d_sh));
	const double s    = x + 0.0066;
	const double env  = a0 + s * (a1 + s * (a2 + s * a3)) + A_bl * std::exp(x / d_bl);
	return nu_inf + step * std::max(0.0, env);
}

// Evaluates RAM-C II 61 km stagnation-line nu_c(x) [1/s] (Liu closure).
inline double ramc61kmStagnationNuC_Liu(const mfem::Vector& p)
{
	const double x = p(0);
	if (x > 0.0) {
		return 0.0;
	}
	constexpr double nu_inf =  3.349800e+07;
	constexpr double x_sh   = -7.241507e-03;
	constexpr double d_sh   =  5.581643e-04;
	constexpr double a0     =  9.148027e+09;
	constexpr double a1     =  1.061459e+12;
	constexpr double a2     = -2.426263e+14;
	constexpr double a3     =  1.968110e+16;
	constexpr double d_bl   =  3.390816e-04;
	constexpr double A_bl   =  1.492221e+09;

	const double step = 0.5 * (1.0 + std::tanh((x - x_sh) / d_sh));
	const double s    = x + 0.0066;
	const double env  = a0 + s * (a1 + s * (a2 + s * a3)) - A_bl * std::exp(x / d_bl);
	return nu_inf + step * std::max(0.0, env);
}

enum class Ramc61kmPreset {
	OmegaP,
	NuC_Gupta,
	NuC_Liu
};

class Ramc61kmStagnationProfile : public Function {
public:
	static constexpr double kNativeLength = 0.012;

	Ramc61kmStagnationProfile(
		Ramc61kmPreset preset,
		double scale = 1.0,
		double nose_x = 0.0,
		int axis = 0,
		double direction = 1.0)
		: preset_(preset),
		  scale_(scale),
		  nose_x_(nose_x),
		  axis_(axis),
		  direction_(direction)
	{
	}

	int dimension() const override { return axis_ + 1; }

	std::unique_ptr<Function> clone() const override
	{
		return std::make_unique<Ramc61kmStagnationProfile>(*this);
	}

	double eval(const Position& pos) const override
	{
		assert(axis_ >= 0 && axis_ < pos.Size());
		mfem::Vector shifted(1);
		shifted(0) = direction_ * (pos[axis_] - nose_x_);
		return evalNative(shifted);
	}

	double evalAtUnitCoord(double xi) const override
	{
		const double u = (direction_ >= 0.0) ? xi : (1.0 - xi);
		mfem::Vector s(1);
		s(0) = -kNativeLength * (1.0 - u);
		return evalNative(s);
	}

private:
	double evalNative(const mfem::Vector& s) const
	{
		switch (preset_) {
		case Ramc61kmPreset::OmegaP:
			return scale_ * ramc61kmStagnationOmegaP(s);
		case Ramc61kmPreset::NuC_Gupta:
			return scale_ * ramc61kmStagnationNuC_Gupta(s);
		case Ramc61kmPreset::NuC_Liu:
			return scale_ * ramc61kmStagnationNuC_Liu(s);
		}
		return 0.0;
	}

	Ramc61kmPreset preset_;
	double scale_;
	double nose_x_;
	int axis_;
	double direction_;
};

// RAM-C II 50 km, Mach 9.5, 12:55 stagnation line.
// x in metres, nose at x = 0, inflow at x = -0.024 m.
// Relative RMS against the H1 axis nodes: 3.6% (largest miss 11% at x = -7.1 mm).
inline double ramc50kmStagnationOmegaP(const mfem::Vector& p)
{
	const double x = p(0);
	if (x > 0.0) {
		return 0.0;
	}
	constexpr double wp_inf =  1.891488933e+05;
	constexpr double x_sh   = -1.404091339e-02;
	constexpr double d_sh   =  3.346945412e-03;
	constexpr double a0     = -1.442830910e+08;
	constexpr double a1     = -8.139327178e+09;
	constexpr double a2     =  2.023092321e+12;
	constexpr double a3     =  1.607015069e+14;
	constexpr double d_bl   =  3.531895113e-03;
	constexpr double A_bl   =  5.860432868e+09;

	const double step = 0.5 * (1.0 + std::tanh((x - x_sh) / d_sh));
	const double s    = x + 1.018700058e-02;
	const double env  = a0 + s * (a1 + s * (a2 + s * a3)) + A_bl * std::exp(x / d_bl);
	return wp_inf + step * std::max(0.0, env);
}

// Relative RMS 1.1% (largest miss 5.9% at x = -20.6 mm).
inline double ramc50kmStagnationNuC_Gupta(const mfem::Vector& p)
{
	const double x = p(0);
	if (x > 0.0) {
		return 0.0;
	}
	constexpr double nu_inf =  1.007955905e+08;
	constexpr double x_sh   = -1.878338741e-02;
	constexpr double d_sh   =  1.433905287e-03;
	constexpr double a0     =  5.432719481e+09;
	constexpr double a1     =  3.957922427e+10;
	constexpr double a2     = -4.318630790e+12;
	constexpr double a3     =  3.187292185e+14;
	constexpr double d_bl   =  4.727815008e-04;
	constexpr double A_bl   =  2.680914883e+09;

	const double step = 0.5 * (1.0 + std::tanh((x - x_sh) / d_sh));
	const double s    = x - 6.151333725e-04;
	const double env  = a0 + s * (a1 + s * (a2 + s * a3)) + A_bl * std::exp(x / d_bl);
	return nu_inf + step * std::max(0.0, env);
}

// Relative RMS 1.4% (largest miss 6.7% at x = -20.6 mm).
inline double ramc50kmStagnationNuC_Liu(const mfem::Vector& p)
{
	const double x = p(0);
	if (x > 0.0) {
		return 0.0;
	}
	constexpr double nu_inf =  1.285807127e+08;
	constexpr double x_sh   = -1.849299301e-02;
	constexpr double d_sh   =  9.186956772e-04;
	constexpr double a0     =  3.007385768e+09;
	constexpr double a1     =  1.495225548e+11;
	constexpr double a2     = -1.230767693e+13;
	constexpr double a3     =  4.136508255e+14;
	constexpr double d_bl   =  4.391254946e-04;
	constexpr double A_bl   =  2.531663926e+09;

	const double step = 0.5 * (1.0 + std::tanh((x - x_sh) / d_sh));
	const double s    = x + 1.698257593e-02;
	const double env  = a0 + s * (a1 + s * (a2 + s * a3)) + A_bl * std::exp(x / d_bl);
	return nu_inf + step * std::max(0.0, env);
}

enum class Ramc50kmPreset {
	OmegaP,
	NuC_Gupta,
	NuC_Liu
};

class Ramc50kmStagnationProfile : public Function {
public:
	static constexpr double kNativeLength = 0.024;

	Ramc50kmStagnationProfile(
		Ramc50kmPreset preset,
		double scale = 1.0,
		double nose_x = 0.0,
		int axis = 0,
		double direction = 1.0)
		: preset_(preset),
		  scale_(scale),
		  nose_x_(nose_x),
		  axis_(axis),
		  direction_(direction)
	{
	}

	int dimension() const override { return axis_ + 1; }

	std::unique_ptr<Function> clone() const override
	{
		return std::make_unique<Ramc50kmStagnationProfile>(*this);
	}

	double eval(const Position& pos) const override
	{
		assert(axis_ >= 0 && axis_ < pos.Size());
		mfem::Vector shifted(1);
		shifted(0) = direction_ * (pos[axis_] - nose_x_);
		return evalNative(shifted);
	}

	double evalAtUnitCoord(double xi) const override
	{
		const double u = (direction_ >= 0.0) ? xi : (1.0 - xi);
		mfem::Vector s(1);
		s(0) = -kNativeLength * (1.0 - u);
		return evalNative(s);
	}

private:
	double evalNative(const mfem::Vector& s) const
	{
		switch (preset_) {
		case Ramc50kmPreset::OmegaP:
			return scale_ * ramc50kmStagnationOmegaP(s);
		case Ramc50kmPreset::NuC_Gupta:
			return scale_ * ramc50kmStagnationNuC_Gupta(s);
		case Ramc50kmPreset::NuC_Liu:
			return scale_ * ramc50kmStagnationNuC_Liu(s);
		}
		return 0.0;
	}

	Ramc50kmPreset preset_;
	double scale_;
	double nose_x_;
	int axis_;
	double direction_;
};

/// Spread σ of exp(-(t-t0)^2 / (2 σ^2)) whose incident *power* is 1/e at `f_1e_hz`.
/// Time is light-metres: |A(ω)|² / |A(0)|² = exp(-(ω σ)²) with ω = 2 π f / c_SI.
inline double gaussianSpreadForPower1e(double f_1e_hz)
{
	if (!(f_1e_hz > 0.0) || !std::isfinite(f_1e_hz)) {
		throw std::runtime_error("f_1e must be > 0 (Hz).");
	}
	return physicalConstants::speedOfLight_SI / (2.0 * M_PI * f_1e_hz);
}

/// Spread σ of exp(-(t-t0)^2 / (2 σ^2)) whose spectrum is `db_cut` dB at `f_max`.
/// |ĝ(f)/ĝ(0)| = exp(-2 π^2 σ^2 f^2) = 10^(db_cut/20), with db_cut < 0.
inline double gaussianSpreadForDbCut(double f_max, double db_cut)
{
	if (!(f_max > 0.0) || !std::isfinite(f_max)) {
		throw std::runtime_error("gaussianSpreadForDbCut: f_max must be > 0.");
	}
	if (!(db_cut < 0.0) || !std::isfinite(db_cut)) {
		throw std::runtime_error("gaussianSpreadForDbCut: db_cut must be < 0.");
	}
	const double pi = std::acos(-1.0);
	const double ratio = std::pow(10.0, db_cut / 20.0);
	const double neg_ln_ratio = -std::log(ratio);
	return std::sqrt(neg_ln_ratio / (2.0 * pi * pi)) / f_max;
}

/// Time signal used as a delta-gap flux. `derivative == false` is the Gaussian
/// pulse. `derivative == true` is its time derivative, so the accumulated field
/// follows the pulse and returns to zero after it passes.
inline double gaussianTimeSignal(
	double t, double t0, double spread, double magnitude, bool derivative)
{
	const double arg = (t - t0) / spread;
	const double pulse = std::exp(-0.5 * arg * arg);
	if (!derivative) {
		return magnitude * pulse;
	}
	return magnitude * (-arg / spread) * pulse;
}

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

	double spread() const { return spread_; }
	const Position& mean() const { return mean_; }
	double frequency() const { return freq_; }

private:
	double spread_;
	Position mean_;
	double freq_;
	int dimension_;
};

/// Incident waveform read from a two-column file: time in seconds, amplitude.
/// Linear interpolation. Zero outside the first and last sample.
/// `eval` receives the same argument as `Gaussian`: phase minus solver time, in light-metres.
class SampledExcitation : public Function {
public:
	explicit SampledExcitation(const std::string& path)
	{
		std::ifstream in(path);
		if (!in) {
			throw std::runtime_error("Cannot open excitation file: " + path);
		}
		std::string line;
		while (std::getline(in, line)) {
			if (line.empty() || line[0] == '#') {
				continue;
			}
			std::istringstream row(line);
			double time = 0.0;
			double amplitude = 0.0;
			if (!(row >> time >> amplitude)) {
				continue;
			}
			if (!std::isfinite(time) || !std::isfinite(amplitude)) {
				throw std::runtime_error("Excitation file has a non-finite sample: " + path);
			}
			if (!times_.empty() && !(time > times_.back())) {
				throw std::runtime_error("Excitation file times must increase: " + path);
			}
			times_.push_back(time);
			values_.push_back(amplitude);
		}
		if (times_.size() < 2) {
			throw std::runtime_error("Excitation file needs at least two samples: " + path);
		}
	}

	int dimension() const { return 1; }

	std::unique_ptr<Function> clone() const {
		return std::make_unique<SampledExcitation>(*this);
	}

	double eval(const Position& pos) const {
		const double time_s = -pos[0] / physicalConstants::speedOfLight_SI;
		if (time_s < times_.front() || time_s > times_.back()) {
			return 0.0;
		}
		const auto upper = std::upper_bound(times_.begin(), times_.end(), time_s);
		if (upper == times_.begin()) {
			return values_.front();
		}
		if (upper == times_.end()) {
			return values_.back();
		}
		const std::size_t i1 = static_cast<std::size_t>(upper - times_.begin());
		const std::size_t i0 = i1 - 1;
		const double span = times_[i1] - times_[i0];
		const double weight = (time_s - times_[i0]) / span;
		return values_[i0] + weight * (values_[i1] - values_[i0]);
	}

private:
	std::vector<double> times_;
	std::vector<double> values_;
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
	/// amplitude_peak: desired max |E| (equatorial |E_θ|) at peak_radius.
	/// Internal scale is set so the analytic Hertzian pulse hits that peak.
	DerivGaussDipole(const double length, const double gaussianSpread, const double gaussMean,
	                 const double amplitude_peak = 1.0, const double peak_radius = 1.0) :
		len_(length),
		gaussSpread_(gaussianSpread),
		gaussMean_(gaussMean),
		amplitude_peak_(amplitude_peak),
		peak_radius_(peak_radius)
	{
		const double unit_peak = peakAbsEtWithUnitScale(
			len_, gaussSpread_, gaussMean_, peak_radius_);
		if (unit_peak <= 0.0 || !(unit_peak == unit_peak)) {
			throw std::runtime_error(
				"DerivGaussDipole: could not compute a positive unit peak |E| "
				"at peak_radius for amplitude_peak normalization.");
		}
		field_scale_ = amplitude_peak_ / unit_peak;
	}

	DerivGaussDipole(const DerivGaussDipole& rhs) :
		len_{ rhs.len_ },
		gaussSpread_{ rhs.gaussSpread_ },
		gaussMean_{ rhs.gaussMean_ },
		amplitude_peak_{ rhs.amplitude_peak_ },
		peak_radius_{ rhs.peak_radius_ },
		field_scale_{ rhs.field_scale_ }
	{
	}

	std::unique_ptr<EHFieldFunction> clone() const {
		return std::make_unique<DerivGaussDipole>(*this);
	}

	double spread() const { return gaussSpread_; }

	double eval(
		const Position& p, const Time& t,
		const FieldType& ft, const Direction& d) const {

		SphericalVector pos(p);

		const auto& cs = physicalConstants::speedOfLight;
		auto cs2 = cs * cs;

		auto radius = pos.radius / cs;
		auto radius2 = radius * radius;
		auto radius3 = radius2 * radius;

		auto sint = std::sin(pos.theta);
		auto cost = std::cos(pos.theta);

		auto spreadsqrt2 = gaussSpread_ * std::sqrt(2.0) ;

		auto expArg = (t - gaussMean_ - pos.radius / cs) / spreadsqrt2;
		auto expArg2 = expArg * expArg;

		auto scalingFactor = (spreadsqrt2 * std::exp(1.0) / 2.0) * field_scale_;

		auto iret = scalingFactor * std::exp(-expArg2);
		auto diret = -iret * 2.0 * expArg / spreadsqrt2;
		auto doublediret = diret * (-2.0) * expArg / spreadsqrt2 + iret * (-2.0) / spreadsqrt2 / spreadsqrt2;

		const auto& ifpe0 = physicalConstants::invFourPiEps0;
		const auto& ifp = physicalConstants::invFourPi;

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

	double length() const { return len_; }
	double gaussianSpread() const { return gaussSpread_; }
	double gaussianMean() const { return gaussMean_; }
	double amplitudePeak() const { return amplitude_peak_; }
	double peakRadius() const { return peak_radius_; }

private:
	/// Max |E_θ| on the equator (θ=π/2) over the pulse, with field_scale = 1.
	static double peakAbsEtWithUnitScale(
		double length, double spread, double mean, double peak_radius)
	{
		using namespace physicalConstants;
		const double cs = speedOfLight;
		const double cs2 = cs * cs;
		const double radius = peak_radius / cs;
		const double radius2 = radius * radius;
		const double radius3 = radius2 * radius;
		const double spreadsqrt2 = spread * std::sqrt(2.0);
		const double scalingFactor = (spreadsqrt2 * std::exp(1.0) / 2.0);
		const double ifpe0 = invFourPiEps0;
		const double sint = 1.0;

		const double t_center = mean + peak_radius / cs;
		const double t_half = 6.0 * spreadsqrt2;
		const int n_samples = 401;
		double peak = 0.0;
		for (int i = 0; i < n_samples; ++i) {
			const double t = t_center - t_half
				+ (2.0 * t_half) * static_cast<double>(i) / static_cast<double>(n_samples - 1);
			const double expArg = (t - mean - peak_radius / cs) / spreadsqrt2;
			const double expArg2 = expArg * expArg;
			const double iret = scalingFactor * std::exp(-expArg2);
			const double diret = -iret * 2.0 * expArg / spreadsqrt2;
			const double doublediret =
				diret * (-2.0) * expArg / spreadsqrt2
				+ iret * (-2.0) / spreadsqrt2 / spreadsqrt2;
			const double et = (length / cs) * ifpe0 * sint
				* (iret / radius3 + diret / (cs * radius2) + doublediret / (cs2 * radius));
			peak = std::max(peak, std::abs(et));
		}
		return peak;
	}

	double len_;
	double gaussSpread_;
	double gaussMean_;
	double amplitude_peak_;
	double peak_radius_;
	double field_scale_; ///< amplitude_peak_ / unit peak |E_θ| at peak_radius_

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

	const Function* function() const { return function_.get(); }
	const Polarization& polarization() const { return polarization_; }
	const Propagation& propagation() const { return propagation_; }
	FieldType fieldType() const { return fieldtype_; }

private:
	std::unique_ptr<Function> function_;
	Polarization polarization_;
	Propagation propagation_;
	FieldType fieldtype_;
};

/// Circular-coax TEM pair. `magnitude` is the conductor voltage \(V_0\), not a uniform \(|E|\).
/// \(\hat{s}\) points from the scattered-field side into the total-field side.
/// \(E_\rho = V_0/(\rho\ln(b/a))\) times the time signal, and \(\mathbf{H}=\hat{s}\times\mathbf{E}\).
class CoaxialMode : public EHFieldFunction {
public:
	CoaxialMode(
		double magnitude,
		double spread,
		double t0,
		bool derivative,
		const mfem::Vector& center,
		const mfem::Vector& axis,
		double inner_radius,
		double outer_radius,
		int total_field_volume,
		int scattered_field_volume) :
		magnitude_{ magnitude },
		spread_{ spread },
		t0_{ t0 },
		derivative_{ derivative },
		inner_radius_{ inner_radius },
		outer_radius_{ outer_radius },
		total_field_volume_{ total_field_volume },
		scattered_field_volume_{ scattered_field_volume }
	{
		if (!(magnitude_ > 0.0) || !std::isfinite(magnitude_)) {
			throw std::runtime_error("coaxial_port magnitude must be > 0.");
		}
		if (!(spread_ > 0.0) || !std::isfinite(spread_)) {
			throw std::runtime_error("coaxial_port spread must be > 0.");
		}
		if (!std::isfinite(t0_)) {
			throw std::runtime_error("coaxial_port t0 must be finite.");
		}
		if (!(outer_radius_ > inner_radius_) || !(inner_radius_ > 0.0)) {
			throw std::runtime_error("coaxial_port radii must satisfy b > a > 0.");
		}
		if (center.Size() != 3 || axis.Size() != 3) {
			throw std::runtime_error("coaxial_port center and axis must be 3-vectors.");
		}
		center_.SetSize(3);
		axis_.SetSize(3);
		for (int d = 0; d < 3; ++d) {
			center_[d] = center[d];
			axis_[d] = axis[d];
		}
		const double n = axis_.Norml2();
		if (!(n > 0.0)) {
			throw std::runtime_error("coaxial_port axis must be nonzero.");
		}
		axis_ /= n;
		log_b_over_a_ = std::log(outer_radius_ / inner_radius_);
	}

	std::unique_ptr<EHFieldFunction> clone() const override {
		return std::make_unique<CoaxialMode>(*this);
	}

	int totalFieldVolume() const { return total_field_volume_; }
	int scatteredFieldVolume() const { return scattered_field_volume_; }
	double innerRadius() const { return inner_radius_; }
	double outerRadius() const { return outer_radius_; }
	double magnitude() const { return magnitude_; }
	double spread() const { return spread_; }

	double eval(
		const Position& p, const Time& t,
		const FieldType& ft, const Direction& d) const override
	{
		double rel[3] = {0.0, 0.0, 0.0};
		for (int i = 0; i < p.Size() && i < 3; ++i) {
			rel[i] = p[i] - center_[i];
		}
		double along = 0.0;
		for (int i = 0; i < 3; ++i) {
			along += rel[i] * axis_[i];
		}
		double radial[3];
		double rho2 = 0.0;
		for (int i = 0; i < 3; ++i) {
			radial[i] = rel[i] - along * axis_[i];
			rho2 += radial[i] * radial[i];
		}
		const double rho = std::sqrt(rho2);
		if (!(rho > 0.0)) {
			return 0.0;
		}
		const double e_rho = gaussianTimeSignal(t, t0_, spread_, magnitude_, derivative_)
			/ (rho * log_b_over_a_);
		double e[3];
		for (int i = 0; i < 3; ++i) {
			e[i] = e_rho * radial[i] / rho;
		}
		if (ft == E) {
			return e[d];
		}
		const double h[3] = {
			axis_[1] * e[2] - axis_[2] * e[1],
			axis_[2] * e[0] - axis_[0] * e[2],
			axis_[0] * e[1] - axis_[1] * e[0]
		};
		return h[d];
	}

private:
	double magnitude_;
	double spread_;
	double t0_;
	bool derivative_;
	double inner_radius_;
	double outer_radius_;
	double log_b_over_a_ = 0.0;
	int total_field_volume_;
	int scattered_field_volume_;
	mfem::Vector center_;
	mfem::Vector axis_;
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
