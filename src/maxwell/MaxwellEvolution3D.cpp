#include "MaxwellEvolution3D.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution3D::MaxwellEvolution3D(
	FiniteElementSpace& fes, Model& model, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{model},
	opts_{ options }
{
	for (auto d: {X, Y, Z}) {
		for (auto f : {E, H}) {
			const auto f2{ altField(f) };
			MS_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(d, fes_), fes_);
			MF_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, d, false, model_, fes_, opts_), fes_);
			MP_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, d, true, model_, fes_, opts_), fes_);
		}
	}
}

void MaxwellEvolution3D::Mult(const Vector& in, Vector& out) const
{
	std::array<Vector, 3> eOld, hOld;
	std::array<GridFunction, 3> eNew, hNew;
	for (int d = X; d <= Z; d++) {
		eOld[d].SetDataAndSize(in.GetData() + d * fes_.GetNDofs(), fes_.GetNDofs());
		hOld[d].SetDataAndSize(in.GetData() + (d + 3) * fes_.GetNDofs(), fes_.GetNDofs());
		eNew[d].MakeRef(&fes_, &out[d * fes_.GetNDofs()]);
		hNew[d].MakeRef(&fes_, &out[(d + 3) * fes_.GetNDofs()]);
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		// dtE_x = MS_y * H_z - MF_y * {H_z} - MP_E * [E_z] +
		//        -MS_z * H_y + MF_z * {H_y} + MP_E * [E_y]
		// Update E.
		MS_[E][z].Mult   (hOld[y], eNew[x]);
		MF_[E][z].AddMult(hOld[y], eNew[x], -1.0);
		MP_[E][z].AddMult(eOld[y], eNew[x], -1.0);
		MS_[E][y].AddMult(hOld[z], eNew[x], -1.0);
		MF_[E][y].AddMult(hOld[z], eNew[x],  1.0);
		MP_[E][y].AddMult(eOld[z], eNew[x],  1.0); 

		// Update H.
		MS_[H][y].Mult   (eOld[z], hNew[x]);
		MF_[H][y].AddMult(eOld[z], hNew[x], -1.0);
		MP_[H][y].AddMult(hOld[z], hNew[x], -1.0);
		MS_[H][z].AddMult(eOld[y], hNew[x], -1.0);
		MF_[H][z].AddMult(eOld[y], hNew[x],  1.0);
		MP_[H][z].AddMult(hOld[y], hNew[x],  1.0);
	}

}

}

