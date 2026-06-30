#pragma once

namespace maxwell {

/// Runtime sign experiments for PML assembly (env PML_SIGN_TEST).
/// 0=default, 1=S1..6 per docs/pml plan sign audit.
enum class PMLSignTestMode {
	Default = 0,
	FlipCorrections = 1,      // S1: swap H/E correction placement signs
	FaceSameAsVolume = 2,     // S2: face driver +w (no SBP opposition)
	FaceCrossColumn = 3,      // S3: face driver uses global one-normal column layout
	NegateDriverWeight = 4,   // S4: negate curl driver weight w
	NegatePMLOperatorMult = 5,// S5: AddMult(..., -1)
	FlipPsiMassSign = 6,      // S6: psi mass placement +1 instead of -1
	IncludeOuterBdyFace = 7     // S7/S8: add terminating boundary face on ψ driver (experimental)
};

PMLSignTestMode getPMLSignTestMode();
const char* pmlSignTestModeName(PMLSignTestMode mode);
double getPMLOperatorMultSign();

} // namespace maxwell
