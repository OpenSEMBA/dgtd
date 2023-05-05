#include "Evolution2D_Spectral.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution2D_Spectral::MaxwellEvolution2D_Spectral(
	FiniteElementSpace& fes, Model& model, SourcesManager& srcmngr, MaxwellEvolOptions& options) :
	TimeDependentOperator(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs()),
	fes_{ fes },
	model_{ model },
	srcmngr_{ srcmngr },
	opts_{ options }
{

	global_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());

	forcing_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());

	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(Y, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { X,Z }, -1.0); //MS
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(X, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { Y,Z });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(X, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { Z,Y });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(Y, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { Z,X }, -1.0);

	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { Y }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { X,Z }); //MFN
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { X }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { Y,Z }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(H, { X }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { Z,Y }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E,	model_, fes_), *buildFluxOperator(H, { Y }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { Z,X });

	if (opts_.fluxType == FluxType::Upwind) {

		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyOperator(H, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { Y }, -1.0); //MP
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyOperator(H, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { X }, -1.0);
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildPenaltyOperator(E, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { E,E }, { Z }, -1.0);

		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { X, X }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { X,X }); //MPNN
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { X, Y }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { X,Y });
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { Y, X }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { Y,X });
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { Y, Y }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { Y,Y });
	
	}

	if (opts_.eigenvals == true) {
		calculateEigenvalues(global_, eigenvals_);
		findMaxEigenvalueModulus(eigenvals_);
	}

	if (opts_.powerMethod != 0)
	{
		pmEigenvalue_ = usePowerMethod(global_, opts_.powerMethod);
	}

	if (opts_.marketFile == true) {
		exportSparseToMarketFile(global_);
	}
}

void MaxwellEvolution2D_Spectral::Mult(const Vector& in, Vector& out) const
{
	out = toMFEMVector( global_ * toEigenVector(in));
}

}

