#include "MaxwellEvolution3D.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution3D::MaxwellEvolution3D(
	FiniteElementSpace& fes, Model& model, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	opts_{ options }
{
	for (auto f : { E, H }) {
		MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyOperator(f, {}, model_, fes_, opts_), fes_);
		for (auto d : { X, Y, Z }) {
			MS_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(d, fes_), fes_);
			for (auto d2 : { X,Y,Z }) {
				for (auto f2 : { E, H }) {
					MFN_[f][f2][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, {d}, model_, fes_), fes_);
					MFNN_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, {d, d2}, model_, fes_), fes_);
				}
			}
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
		eNew[d] = 0.0;
		hNew[d] = 0.0;
	}

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;

		//Centered
		MS_[H][y]		 ->AddMult(eOld[z], hNew[x], -1.0);
		MS_[H][z]		 ->AddMult(eOld[y], hNew[x]);
		MFN_[H][E][y]  	 ->AddMult(eOld[z], hNew[x], 1.0);
		MFN_[H][E][z]	 ->AddMult(eOld[y], hNew[x], -1.0);
		
		MS_[E][y]		 ->AddMult(hOld[z], eNew[x]);
		MS_[E][z]		 ->AddMult(hOld[y], eNew[x], -1.0);
		MFN_[E][H][y]  	 ->AddMult(hOld[z], eNew[x], -1.0);
		MFN_[E][H][z]	 ->AddMult(hOld[y], eNew[x], 1.0);

		if (opts_.fluxType == FluxType::Upwind) {
			MFNN_[H][H][X][x]->AddMult(hOld[X], hNew[x], 1.0);
			MFNN_[H][H][Y][x]->AddMult(hOld[Y], hNew[x], 1.0);
			MFNN_[H][H][Z][x]->AddMult(hOld[Z], hNew[x], 1.0);
			MP_[H]			 ->AddMult(hOld[x], hNew[x],-1.0);
			
			MFNN_[E][E][X][x]->AddMult(eOld[X], eNew[x], 1.0);
			MFNN_[E][E][Y][x]->AddMult(eOld[Y], eNew[x], 1.0);
			MFNN_[E][E][Z][x]->AddMult(eOld[Z], eNew[x], 1.0);
			MP_[E]			 ->AddMult(eOld[x], eNew[x],-1.0);
		}

	}


}

}

