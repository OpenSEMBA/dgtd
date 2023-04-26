#include "Evolution1D_Spectral.h"


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

	global_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
				   numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());

	allocateDenseInEigen1D(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(X, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, -1.0); // MS
	allocateDenseInEigen1D(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(X, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, -1.0);

	allocateDenseInEigen1D(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { X }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }); // MF
	allocateDenseInEigen1D(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(H, { X }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H });

	if (opts_.fluxType == FluxType::Upwind) {
		allocateDenseInEigen1D(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyOperator(H, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, -1.0); //MP
		allocateDenseInEigen1D(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildPenaltyOperator(E, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { E,E }, -1.0);
	}

	forcing_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<Planewave*>(source.get())) {
			allocateDenseInEigen1D(buildIBFIByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxFunctionOperator1D(model_, fes_), fes_)->SpMat().ToDenseMatrix(), forcing_, { H,H }); //MBF
			allocateDenseInEigen1D(buildIBFIByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxFunctionOperator1D(model_, fes_), fes_)->SpMat().ToDenseMatrix(), forcing_, { E,E });
			if (opts_.fluxType == FluxType::Upwind) {
				allocateDenseInEigen1D(buildIBFIByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyFunctionOperator1D(model_, fes_), fes_)->SpMat().ToDenseMatrix(), forcing_, { H,H }); //MBP
				allocateDenseInEigen1D(buildIBFIByMult(*buildInverseMassMatrix(E, model_, fes_), *buildPenaltyFunctionOperator1D(model_, fes_), fes_)->SpMat().ToDenseMatrix(), forcing_, { E,E });
			}
		}
	}

	if (opts_.eigenvals == true) {
		SparseMatrix mat{ toMFEMSparse(global_) };
		calculateEigenvalues(mat, eigenvals_);
		//checkEigenvalues(eigenvals_);
		eigenvals_.Print(std::cout);
	}

	if (opts_.powerMethod != 0)
	{
		pmEigenvalue_ = usePowerMethod(global_, opts_.powerMethod);
		std::cout << "And the powerMethoded Eigenval is..." + std::to_string(pmEigenvalue_) << std::endl;
	}

	if (opts_.marketFile == true) {
		exportSparseToMarketFile(global_);
	}

}


void MaxwellEvolution1D_Spectral::Mult(const Vector& in, Vector& out) const
{
	Eigen::VectorXd fieldsOld{ toEigenVector(in) };
	Eigen::VectorXd fieldsNew{ global_ * fieldsOld };

	for (const auto& source : srcmngr_.sources) {
		if (dynamic_cast<Planewave*>(source.get())) {
			std::array<std::array<GridFunction, 3>, 2> eFunc(srcmngr_.evalTimeVarField(GetTime()));
			std::array<std::array<GridFunction, 3>, 2> hFunc(srcmngr_.evalTimeVarField(GetTime()));

			Eigen::VectorXd forcVecsOld;
			forcVecsOld.resize(2 * eFunc[E][Y].Size());
			for (int i = 0; i < eFunc[E][Y].Size(); ++i) {
				forcVecsOld(i)                      = eFunc[E][Y].Elem(i);
				forcVecsOld(i + eFunc[E][Y].Size()) = hFunc[H][Z].Elem(i);
			} 
			fieldsNew += forcing_ * forcVecsOld;
		}
	}
	out = toMFEMVector(fieldsNew);
}

}

