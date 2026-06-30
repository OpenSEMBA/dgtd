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
	switch (kind_) {
	case Kind::Alpha:
		return prof.alpha;
	case Kind::Sigma:
		return prof.sigma;
	case Kind::InvKappa:
		return (prof.kappa > 0.0) ? (1.0 / prof.kappa) : 0.0;
	}
	return 0.0;
}

} // namespace maxwell
