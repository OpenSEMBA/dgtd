#include "PMLDGHelpers.h"

namespace maxwell {

double pmlCurlDerivWeight(FieldType out_f, Direction out_c, FieldType in_f,
                          Direction in_c, Direction deriv)
{
	if (out_f == E && in_f == H) {
		if (out_c == X) {
			if (in_c == Z && deriv == Y) return 1.0;
			if (in_c == Y && deriv == Z) return -1.0;
		} else if (out_c == Y) {
			if (in_c == X && deriv == Z) return 1.0;
			if (in_c == Z && deriv == X) return -1.0;
		} else if (out_c == Z) {
			if (in_c == Y && deriv == X) return 1.0;
			if (in_c == X && deriv == Y) return -1.0;
		}
	} else if (out_f == H && in_f == E) {
		if (out_c == X) {
			if (in_c == Y && deriv == Z) return -1.0;
			if (in_c == Z && deriv == Y) return 1.0;
		} else if (out_c == Y) {
			if (in_c == Z && deriv == X) return -1.0;
			if (in_c == X && deriv == Z) return 1.0;
		} else if (out_c == Z) {
			if (in_c == X && deriv == Y) return -1.0;
			if (in_c == Y && deriv == X) return 1.0;
		}
	}
	return 0.0;
}

bool pmlPsiEComponentActive(Direction h_out_c, Direction stretch_dir)
{
	for (Direction in_c = X; in_c <= Z; ++in_c) {
		if (pmlCurlDerivWeight(H, h_out_c, E, in_c, stretch_dir) != 0.0) {
			return true;
		}
	}
	return false;
}

bool pmlPsiHComponentActive(Direction e_out_c, Direction stretch_dir)
{
	for (Direction in_c = X; in_c <= Z; ++in_c) {
		if (pmlCurlDerivWeight(E, e_out_c, H, in_c, stretch_dir) != 0.0) {
			return true;
		}
	}
	return false;
}

bool pmlPsiEDriverCoupling(Direction h_out_c, Direction stretch_dir,
                           Direction& in_c, double& weight)
{
	for (Direction ec = X; ec <= Z; ++ec) {
		const double w = pmlCurlDerivWeight(H, h_out_c, E, ec, stretch_dir);
		if (w != 0.0) {
			in_c = ec;
			weight = w;
			return true;
		}
	}
	return false;
}

bool pmlPsiHDriverCoupling(Direction e_out_c, Direction stretch_dir,
                           Direction& in_c, double& weight)
{
	for (Direction hc = X; hc <= Z; ++hc) {
		const double w = pmlCurlDerivWeight(E, e_out_c, H, hc, stretch_dir);
		if (w != 0.0) {
			in_c = hc;
			weight = w;
			return true;
		}
	}
	return false;
}

} // namespace maxwell
