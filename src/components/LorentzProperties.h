#pragma once

#include "Types.h"
#include "math/Function.h"

#include <memory>
#include <nlohmann/json.hpp>
#include <vector>

namespace maxwell {

/// Single-pole electric Lorentz (or Drude when omega_1 == 0) on one mesh attribute.
/// Rates are SI values divided by c_SI (solver time is light-meters).
/// epsilon_ stored on Material for this tag is eps_inf.
struct LorentzProperties {
	Attribute geom_tag = 0;
	double eps_inf = 1.0;
	double omega_p = 0.0;
	double omega_1 = 0.0;
	double gamma = 0.0;

	std::shared_ptr<const Function> omega_p_profile;
	std::shared_ptr<const Function> omega_1_profile;
	std::shared_ptr<const Function> gamma_profile;

	/// 1D physical bounding interval [region_min_x, region_max_x] of the mesh region
	/// carrying this material. When region_max_x > region_min_x, 1D profiles stretch
	/// across [region_min_x, region_max_x] via normalized coordinate xi in [0, 1].
	double region_min_x = 0.0;
	double region_max_x = 0.0;

	bool hasRegionBounds() const
	{
		return region_max_x > region_min_x;
	}

	bool hasSpatialProfile() const
	{
		return omega_p_profile != nullptr ||
		       omega_1_profile != nullptr ||
		       gamma_profile != nullptr;
	}

	double evalProfile(const Function& fn, const mfem::Vector& x) const
	{
		if (x.Size() == 1 && hasRegionBounds()) {
			const double xi = (x(0) - region_min_x) / (region_max_x - region_min_x);
			return fn.evalAtUnitCoord(xi);
		}
		return fn.eval(x);
	}

	double evalOmegaP(const mfem::Vector& x) const
	{
		return omega_p_profile ? evalProfile(*omega_p_profile, x) : omega_p;
	}

	double evalOmega1(const mfem::Vector& x) const
	{
		return omega_1_profile ? evalProfile(*omega_1_profile, x) : omega_1;
	}

	double evalGamma(const mfem::Vector& x) const
	{
		return gamma_profile ? evalProfile(*gamma_profile, x) : gamma;
	}
};

/// Parse a `lorentz` object. geom_tag is filled by the caller.
LorentzProperties parseLorentzObject(const nlohmann::json& lorentz_json);

/// Parse a `drude` object (cold-plasma Lorentz pole with omega_1 = 0). geom_tag is filled by the caller.
LorentzProperties parseDrudeObject(const nlohmann::json& drude_json);

} // namespace maxwell
