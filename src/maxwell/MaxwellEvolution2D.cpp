#include "MaxwellEvolution2D.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution2D::MaxwellEvolution2D(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{
	for (auto f : { E, H }) {
		if (opts_.fluxType == FluxType::Upwind) {
			MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyOperator(f, {}, model_, fes_, opts_), fes_);
			MBP_[f] = buildIBFIByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyFunctionOperator(f, model_, fes_), fes_);
		}
		for (auto d : { X, Y, Z }) {
			MS_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(d, fes_), fes_);
			for (auto d2 : { X,Y,Z }) {
				for (auto f2 : { E, H }) {
					MFN_[f][f2][d]      = buildByMult    (*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, {d}, model_, fes_), fes_);
					if (opts_.fluxType == FluxType::Upwind) {
						MPNN_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, { d, d2 }, model_, fes_), fes_);
						MBFN_[f][f2][d] = buildIBFIByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxFunctionOperator(f2, { d }, model_, fes_), fes_);
						MBPNN_[f][f][d][d2] = buildIBFIByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyFunctionOperator(f, model_, fes_), fes_);
					}
				}
			}
		}
	}
}

void MaxwellEvolution2D::Mult(const Vector& in, Vector& out) const
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

	MFN_[H][E][Y]    ->AddMult(eOld[Z], hNew[X]);
	MS_[H][Y]        ->AddMult(eOld[Z], hNew[X], -1.0);

	MFN_[H][E][X]    ->AddMult(eOld[Z], hNew[Y], -1.0);				 
	MS_	[H][X]       ->AddMult(eOld[Z], hNew[Y]);

	MFN_[E][H][Y]->AddMult(hOld[X], eNew[Z]);
	MFN_[E][H][X]->AddMult(hOld[Y], eNew[Z], -1.0);

	MS_[E][X]	 ->AddMult(hOld[Y], eNew[Z]);
	MS_[E][Y]    ->AddMult(hOld[X], eNew[Z], -1.0);

	if (opts_.fluxType == FluxType::Upwind) {
		MPNN_[H][H][X][X]->AddMult(hOld[X], hNew[X]);
		MPNN_[H][H][Y][X]->AddMult(hOld[Y], hNew[X]);
		MP_[H]->AddMult(hOld[X], hNew[X], -1.0);

		MPNN_[H][H][X][Y]->AddMult(hOld[X], hNew[Y]);
		MPNN_[H][H][Y][Y]->AddMult(hOld[Y], hNew[Y]);
		MP_[H]->AddMult(hOld[Y], hNew[Y], -1.0);

		MP_[E]->AddMult(eOld[Z], eNew[Z], -1.0);
	}

}

}

