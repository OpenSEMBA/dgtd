#include "MaxwellEvolution1D.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution1D::MaxwellEvolution1D(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{
	for (auto f : {E, H}) {
		const auto f2{ altField(f) };
		MS_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(X, fes_), fes_);
		MF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator1D(f2, {X}, model_, fes_), fes_);
		MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyOperator1D(f, {}, model_, fes_, opts_), fes_);
		MT_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFunctionOperator1D(f, model_, fes_), fes_);
	}
	buildFunctionOperator1D(E, model_, fes_).get()->SpMat().Print(std::cout);
	std::cout << ".-.-.-." << std::endl;
}

void MaxwellEvolution1D::Mult(const Vector& in, Vector& out) const
{
	Vector eOld, hOld;
	GridFunction eNew, hNew;

	eOld.SetDataAndSize(in.GetData(), fes_.GetNDofs());
	hOld.SetDataAndSize(in.GetData() + fes_.GetNDofs(), fes_.GetNDofs());
	eNew.MakeRef(&fes_, &out[0]);
	hNew.MakeRef(&fes_, &out[fes_.GetNDofs()]);


	// dtE = - MS * H + MF * [H] (signs in coeff)
	// Update E.
	MF_[E]->Mult   (hOld, eNew);
	MS_[E]->AddMult(hOld, eNew, -1.0);

	// dtH = - MS * E + MF * [E] (signs in coeff)
	// Update H.
	MF_[H]->Mult   (eOld, hNew);
	MS_[H]->AddMult(eOld, hNew, -1.0);

	if (opts_.fluxType == FluxType::Upwind) {
		MP_[E]->AddMult(eOld, eNew, -1.0);
		MP_[H]->AddMult(hOld, hNew, -1.0);
	}

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<PlaneWave*>(source.get())) {
			GridFunction eFunc(srcmngr_.evalTotalField(GetTime()));
			GridFunction hFunc(srcmngr_.evalTotalField(GetTime()));
			MT_[E]->AddMult(eFunc, eNew);
			MT_[H]->AddMult(hFunc, hNew);
		}
	}


}

}

