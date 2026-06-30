#include "PMLSignTest.h"

#include <cstdlib>
#include <cstring>

namespace maxwell {

PMLSignTestMode getPMLSignTestMode()
{
	const char* env = std::getenv("PML_SIGN_TEST");
	if (!env || env[0] == '\0') {
		return PMLSignTestMode::Default;
	}
	const int v = std::atoi(env);
	if (v >= 0 && v <= 8) {
		return static_cast<PMLSignTestMode>(v);
	}
	return PMLSignTestMode::Default;
}

const char* pmlSignTestModeName(PMLSignTestMode mode)
{
	switch (mode) {
	case PMLSignTestMode::Default: return "default";
	case PMLSignTestMode::FlipCorrections: return "S1_flip_corrections";
	case PMLSignTestMode::FaceSameAsVolume: return "S2_face_same_as_vol";
	case PMLSignTestMode::FaceCrossColumn: return "S3_face_cross_column";
	case PMLSignTestMode::NegateDriverWeight: return "S4_negate_driver_w";
	case PMLSignTestMode::NegatePMLOperatorMult: return "S5_negate_pml_mult";
	case PMLSignTestMode::FlipPsiMassSign: return "S6_flip_psi_mass";
	case PMLSignTestMode::IncludeOuterBdyFace: return "S7_include_outer_bdy_face";
	}
	return "unknown";
}

double getPMLOperatorMultSign()
{
	if (getPMLSignTestMode() == PMLSignTestMode::NegatePMLOperatorMult) {
		return -1.0;
	}
	return 1.0;
}

} // namespace maxwell
