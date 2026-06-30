#pragma once

#include "PMLProfiles.h"
#include "Types.h"

namespace maxwell {

/// Evaluate PML profile coefficient at a quadrature point (zero outside PML elements).
class PMLProfileCoefficient : public mfem::Coefficient {
public:
	enum class Kind { Alpha, Sigma, InvKappa };

	PMLProfileCoefficient(const PMLProfileData& profiles, Direction stretch_dir, Kind kind);

	double Eval(mfem::ElementTransformation& T,
	            const mfem::IntegrationPoint& ip) override;

private:
	const PMLProfileData& profiles_;
	Direction stretch_dir_;
	Kind kind_;
};

} // namespace maxwell
