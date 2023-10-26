#include "EigenvalueEstimator.h"

namespace maxwell {

int EigenvalueEstimator::getOffset(const FieldType& ft, const Direction& d)
{
	switch (ft) {
	case E:
		switch (d) {
		case X:
			return 0;
		case Y:
			return fes_.GetNDofs();
		case Z:
			return 2 * fes_.GetNDofs();
		default:
			throw std::runtime_error("Incorrect Direction for FieltType E in getOffset.");
		}
	case H:
		switch (d) {
		case X:
			return 3 * fes_.GetNDofs();
		case Y:
			return 4 * fes_.GetNDofs();
		case Z:
			return 5 * fes_.GetNDofs();
		default:
			throw std::runtime_error("Incorrect Direction for FieltType H in getOffset.");
		}
	default:
		throw std::runtime_error("Incorrect FieldType in getOffset.");
	}
}

EigenvalueEstimator::EigenvalueEstimator(
	FiniteElementSpace& fes,
	Model& model,
	EvolutionOptions& options) :
	fes_{ fes },
	model_{ model },
	opts_{ options }
{

	mat_.resize(numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs(),
		numberOfFieldComponents * numberOfMaxDimensions * fes.GetNDofs());

	mat_.setConstant(0.0);

	auto invM_E{ toEigen(*buildInverseMassMatrix(E, model_, fes_)->SpMat().ToDenseMatrix()) };
	auto invM_H{ toEigen(*buildInverseMassMatrix(H, model_, fes_)->SpMat().ToDenseMatrix()) };
	Eigen::MatrixXd invM(invM_E);

	for (auto f : { E, H }) {

		f == E ? invM = invM_E : invM = invM_H;

		for (int d = X; d < Z; d++) {

			//MP
			mat_.block(getOffset(f, d), getOffset(f, d), fes_.GetNDofs(), fes_.GetNDofs()) +=
				invM * toEigen(*buildZeroNormalOperator(f, model_, fes_, opts_)->SpMat().ToDenseMatrix());
			//MFNN
			mat_.block(getOffset(f, X), getOffset(f, d), fes_.GetNDofs(), fes_.GetNDofs()) -=
				invM * toEigen(*buildTwoNormalOperator(altField(f), { d, X }, model_, fes_, opts_)->SpMat().ToDenseMatrix());
			mat_.block(getOffset(f, Y), getOffset(f, d), fes_.GetNDofs(), fes_.GetNDofs()) -=
				invM * toEigen(*buildTwoNormalOperator(altField(f), { d, Y }, model_, fes_, opts_)->SpMat().ToDenseMatrix());
			mat_.block(getOffset(f, Z), getOffset(f, d), fes_.GetNDofs(), fes_.GetNDofs()) -=
				invM * toEigen(*buildTwoNormalOperator(altField(f), { d, Z }, model_, fes_, opts_)->SpMat().ToDenseMatrix());

		}

		f == H ? invM = -invM_H : invM = invM_E;

		for (int x = X; x <= Z; x++) {
			int y = (x + 1) % 3;
			int z = (x + 2) % 3;

			mat_.block(getOffset(f, x), getOffset(altField(f), y), fes_.GetNDofs(), fes_.GetNDofs()) -=
				invM * toEigen(*buildDerivativeOperator(z, fes_)->SpMat().ToDenseMatrix());
			mat_.block(getOffset(f, x), getOffset(altField(f), y), fes_.GetNDofs(), fes_.GetNDofs()) -=
				invM * toEigen(*buildOneNormalOperator(altField(f), { z }, model_, fes_, opts_)->SpMat().ToDenseMatrix());

			mat_.block(getOffset(f, x), getOffset(altField(f), z), fes_.GetNDofs(), fes_.GetNDofs()) +=
				invM * toEigen(*buildDerivativeOperator(y, fes_)->SpMat().ToDenseMatrix());
			mat_.block(getOffset(f, x), getOffset(altField(f), z), fes_.GetNDofs(), fes_.GetNDofs()) +=
				invM * toEigen(*buildOneNormalOperator(altField(f), { y }, model_, fes_, opts_)->SpMat().ToDenseMatrix());
		}
	}
}
}