#include "PMLCoefficients.h"

#include <cmath>

namespace maxwell {

PMLProfileCoefficient::PMLProfileCoefficient(
	const PMLProfileData& profiles, Direction stretch_dir, Kind kind)
	: profiles_(profiles), stretch_dir_(stretch_dir), kind_(kind)
{
}

double PMLProfileCoefficient::Eval(
	mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip)
{
	PMLDirectionProfiles prof;
	profiles_.evaluateAtTransform(T, ip, stretch_dir_, prof);
	const double inv_kappa = (prof.kappa > 0.0) ? (1.0 / prof.kappa) : 0.0;
	switch (kind_) {
	case Kind::Alpha:
		return prof.alpha;
	case Kind::Sigma:
		return prof.sigma;
	case Kind::InvKappa:
		return inv_kappa;
	case Kind::Decay:
		// Gedney ADE: ∂ψ/∂t = -(α + σ/κ) ψ + (σ/κ) D(F)
		return prof.alpha + prof.sigma * inv_kappa;
	case Kind::SigmaOverKappa:
		return prof.sigma * inv_kappa;
	}
	return 0.0;
}

} // namespace maxwell
