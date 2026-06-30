#pragma once

#include "Types.h"

namespace maxwell {

/// Curl derivative weight for discrete Maxwell (matches DGOperatorFactory signs).
double pmlCurlDerivWeight(FieldType out_f, Direction out_c, FieldType in_f,
                            Direction in_c, Direction deriv);

/// True if ψ^E_{stretch_d, h_out_c} has a curl driver/correction pair.
bool pmlPsiEComponentActive(Direction h_out_c, Direction stretch_dir);

/// True if ψ^H_{stretch_d, e_out_c} has a curl driver/correction pair.
bool pmlPsiHComponentActive(Direction e_out_c, Direction stretch_dir);

/// Input E component and driver weight for ψ^E_{stretch_d, h_out_c}. Returns false if inactive.
bool pmlPsiEDriverCoupling(Direction h_out_c, Direction stretch_dir,
                           Direction& in_c, double& weight);

/// Input H component and driver weight for ψ^H_{stretch_d, e_out_c}. Returns false if inactive.
bool pmlPsiHDriverCoupling(Direction e_out_c, Direction stretch_dir,
                           Direction& in_c, double& weight);

} // namespace maxwell
