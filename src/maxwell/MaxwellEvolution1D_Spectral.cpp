#include "MaxwellEvolution1D_Spectral.h"


namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution1D_Spectral::MaxwellEvolution1D_Spectral(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{
	std::array<FiniteElementOperator, 2> MS, MF, MP;
	std::array<FiniteElementIBFIOperator, 2> MBF, MBP;
	for (auto f : {E, H}) {
		const auto f2{ altField(f) };
		MS[f] = buildByMult		 (*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(X, fes_), fes_);
		MF[f] = buildByMult		 (*buildInverseMassMatrix(f, model_, fes_), *buildFluxOperator1D(f2, {X}, model_, fes_), fes_);
		MBF[f] = buildIBFIByMult	 (*buildInverseMassMatrix(f, model_, fes_), *buildFluxFunctionOperator1D(model_, fes_), fes_);
		if (opts_.fluxType == FluxType::Upwind) {
			MP[f] = buildByMult	 (*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyOperator1D(f, {}, model_, fes_, opts_), fes_);
			MBP[f] = buildIBFIByMult(*buildInverseMassMatrix(f, model_, fes_), *buildPenaltyFunctionOperator1D(model_, fes_), fes_);
		}
	}

	global_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
				   numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());
	global_.setZero();
	allocateDenseInEigen1D<FiniteElementOperator>(MS, global_, -1.0, true);
	allocateDenseInEigen1D<FiniteElementOperator>(MF, global_,  1.0, true);

	if (opts_.fluxType == FluxType::Upwind) {
		allocateDenseInEigen1D<FiniteElementOperator>(MP, global_, -1.0);
	}

	forcing_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());
	forcing_.setZero();

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<PlaneWave*>(source.get())) {
			allocateDenseInEigen1D<FiniteElementIBFIOperator>(MBF, forcing_);
			if (opts_.fluxType == FluxType::Upwind) {
				allocateDenseInEigen1D<FiniteElementIBFIOperator>(MBP, forcing_);
			}
		}
	}

	calculateEigenvalues(global_, eigenvals_);
	checkEigenvalues(eigenvals_);
	exportSparseToMarketFile(global_);

}


void MaxwellEvolution1D_Spectral::Mult(const Vector& in, Vector& out) const
{
	Eigen::VectorXd fieldsOld{ toEigenVector(in) };
	Eigen::VectorXd fieldsNew{ global_ * fieldsOld };

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<PlaneWave*>(source.get())) {
			GridFunction eFunc(srcmngr_.evalTotalField(GetTime()));
			GridFunction hFunc(srcmngr_.evalTotalField(GetTime()));

			Eigen::VectorXd forcVecsOld;
			forcVecsOld.resize(2 * eFunc.Size());
			for (int i = 0; i < eFunc.Size(); ++i) {
				forcVecsOld(i)                = eFunc.Elem(i);
				forcVecsOld(i + eFunc.Size()) = hFunc.Elem(i);
			}
			fieldsNew += forcing_ * forcVecsOld;
		}
	}
	out = toMFEMVector(fieldsNew);
}

}

