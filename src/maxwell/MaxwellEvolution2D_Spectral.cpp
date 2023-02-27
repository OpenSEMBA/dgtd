#include "MaxwellEvolution2D_Spectral.h"

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
	std::array<std::array<						FiniteElementOperator, 3>, 2> MS;
	std::array<std::array<std::array<			FiniteElementOperator, 3>, 2>, 2> MFN;

	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MPNN;
	std::array<									FiniteElementOperator, 2> MP;

	std::array<std::array<std::array<std::array<FiniteElementIBFIOperator, 3>, 3>, 2>, 2> MBPNN;
	std::array<std::array<std::array<			FiniteElementIBFIOperator, 3>, 2>, 2> MBFN;
	std::array<									FiniteElementIBFIOperator, 2> MBP;

	global_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());
	global_.setZero();

	forcing_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());
	forcing_.setZero();

	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(Y, fes_), fes_), global_, { H,E }, { X,Z }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildDerivativeOperator(X, fes_), fes_), global_, { H,E }, { Y,Z });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(X, fes_), fes_), global_, { E,H }, { Z,Y });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildDerivativeOperator(Y, fes_), fes_), global_, { E,H }, { Z,X }, -1.0);

	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { Y }, model_, fes_), fes_), global_, { H,E }, { X,Z });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(E, { X }, model_, fes_), fes_), global_, { H,E }, { Y,Z }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildFluxOperator(H, { X }, model_, fes_), fes_), global_, { E,H }, { Z,Y }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E,	model_, fes_), *buildFluxOperator(H, { Y }, model_, fes_), fes_), global_, { E,H }, { Z,X });

	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyOperator(H, {}, model_, fes_, opts_), fes_), global_, { H,H }, { Y }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildPenaltyOperator(H, {}, model_, fes_, opts_), fes_), global_, { H,H }, { X }, -1.0);
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(E, model_, fes_), *buildPenaltyOperator(E, {}, model_, fes_, opts_), fes_), global_, { E,E }, { Z }, -1.0);

	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { X, X }, model_, fes_), fes_), global_, { H,H }, { X,X });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { X, Y }, model_, fes_), fes_), global_, { H,H }, { X,Y });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { Y, X }, model_, fes_), fes_), global_, { H,H }, { Y,X });
	allocateDenseInEigen(buildByMult(*buildInverseMassMatrix(H, model_, fes_), *buildFluxOperator(H, { Y, Y }, model_, fes_), fes_), global_, { H,H }, { Y,Y });

}

void MaxwellEvolution2D_Spectral::Mult(const Vector& in, Vector& out) const
{
	Eigen::VectorXd fieldsOld{ toEigenVector(in) };
	Eigen::VectorXd fieldsNew{ global_ * fieldsOld };
	out = toMFEMVector(fieldsNew);
}

}

