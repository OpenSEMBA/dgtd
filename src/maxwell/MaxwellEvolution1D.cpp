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
		MF_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, X, false, model_, fes_, opts_), fes_);
		MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator(f2, X, true, model_, fes_, opts_), fes_);
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

	// dtE_x = MS_y * H_z - MF_y * {H_z} - MP_E * [E_z] +
	//        -MS_z * H_y + MF_z * {H_y} + MP_E * [E_y]
	// Update E.
	MS_[E]->Mult   (hOld, eNew);
	MF_[E]->AddMult(hOld, eNew, -1.0);
	MP_[E]->AddMult(eOld, eNew, -1.0);
	MS_[E]->AddMult(hOld, eNew, -1.0);
	MF_[E]->AddMult(hOld, eNew,  1.0);
	MP_[E]->AddMult(eOld, eNew,  1.0); 

	// Update H.
	MS_[H]->Mult   (eOld, hNew);
	MF_[H]->AddMult(eOld, hNew, -1.0);
	MP_[H]->AddMult(hOld, hNew, -1.0);
	MS_[H]->AddMult(eOld, hNew, -1.0);
	MF_[H]->AddMult(eOld, hNew,  1.0);
	MP_[H]->AddMult(hOld, hNew,  1.0);


}

}

