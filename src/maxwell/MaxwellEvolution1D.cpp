#include "MaxwellEvolution1D.h"


namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution1D::MaxwellEvolution1D(
	FiniteElementSpace& fes, Model& model, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	opts_{ options }
{
	for (auto f : {E, H}) {
		const auto f2{ altField(f) };
		MS_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(X, fes_), fes_);
		MF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxJumpOperator(f2, std::vector<Direction>{}, false, model_, fes_, opts_), fes_);
		/*MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxJumpOperator(f2, std::vector<Direction>{}, true, model_, fes_, opts_), fes_);*/
	}
}

void MaxwellEvolution1D::Mult(const Vector& in, Vector& out) const
{
	Vector eOld, hOld;
	GridFunction eNew, hNew;

	eOld.SetDataAndSize(in.GetData(), fes_.GetNDofs());
	hOld.SetDataAndSize(in.GetData() + fes_.GetNDofs(), fes_.GetNDofs());
	eNew.MakeRef(&fes_, &out[0]);
	hNew.MakeRef(&fes_, &out[fes_.GetNDofs()]);


	// dtE = - MS * H + MF * [H] - MF * [E]
	// Update E.
	MF_[E]->Mult   (hOld, eNew);
	MS_[E]->AddMult(hOld, eNew, -1.0);
	MF_[E]->AddMult(eOld, eNew, -1.0);

	// dtH = - MS * E + MF * [E] - MF * [H]
	// Update H.
	MF_[H]->Mult   (eOld, hNew);
	MS_[H]->AddMult(eOld, hNew, -1.0);
	MF_[H]->AddMult(hOld, hNew, -1.0);



}

}

