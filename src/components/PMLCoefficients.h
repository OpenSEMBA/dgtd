#pragma once

#include "PMLProfiles.h"
#include "Types.h"

namespace maxwell {

/// Evaluate PML profile coefficient at a quadrature point (zero outside PML elements).
class PMLProfileCoefficient : public mfem::Coefficient {
public:
	/// Profile coefficient kinds used by ADE CFS-CPML.
	/// Decay = α + σ/κ is the Gedney ADE auxiliary damping (not α alone).
	/// SigmaOverKappa = σ/κ is the consistent driver when field corrections use ψ/κ.
	enum class Kind { Alpha, Sigma, InvKappa, Decay, SigmaOverKappa };

	PMLProfileCoefficient(const PMLProfileData& profiles, Direction stretch_dir, Kind kind);

	double Eval(mfem::ElementTransformation& T,
	            const mfem::IntegrationPoint& ip) override;

private:
	const PMLProfileData& profiles_;
	Direction stretch_dir_;
	Kind kind_;
};

} // namespace maxwell
