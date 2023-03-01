#include "MaxwellEvolution3D_Spectral.h"

namespace maxwell {

using namespace mfem;
using namespace mfemExtension;

MaxwellEvolution3D_Spectral::MaxwellEvolution3D_Spectral(
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

	for (int x = X; x <= Z; x++) {
		int y = (x + 1) % 3;
		int z = (x + 2) % 3;
		
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(y, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { x,z }, -1.0); // MS
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(z, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { x,y });
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(y, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { x,z });
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(z, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { x,y }, -1.0);

		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { y }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { x,z }); // MFN
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { z }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,E }, { x,y }, -1.0);
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(H, { y }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { x,z }, -1.0);
		allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(H, { z }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,H }, { x,y });

		if (opts_.fluxType == FluxType::Upwind) {

			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyOperator(H, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { x }, -1.0); // MP
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildPenaltyOperator(E, {}, model_, fes_, opts_), fes_)->SpMat().ToDenseMatrix(), global_, { E,E }, { x }, -1.0);

			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { X, x }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { X,x }); //MPNN
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { Y, x }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { Y,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { Z, x }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { H,H }, { Z,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(E, { X, x }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,E }, { X,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(E, { Y, x }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,E }, { Y,x });
			allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(E, { Z, x }, model_, fes_), fes_)->SpMat().ToDenseMatrix(), global_, { E,E }, { Z,x });

		}

	}

	if (opts_.eigenvals == true) {
		calculateEigenvalues(global_, eigenvals_);
		checkEigenvalues(eigenvals_);
	}

	if (opts_.marketFile == true) {
		exportSparseToMarketFile(global_);
	}

 }

void MaxwellEvolution3D_Spectral::Mult(const Vector& in, Vector& out) const
{
	out = toMFEMVector(global_ * toEigenVector(in));
}

}

