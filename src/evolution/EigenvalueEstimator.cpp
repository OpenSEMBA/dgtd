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
		throw std::runtime_error("Incorrect Direction in getOffset.");
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

	for (auto f : { E, H }) {
		MP_[f] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildZeroNormalOperator(f, model_, fes_, opts_), fes_);
		for (auto d{ X }; d <= Z; d++) {
			MS_[f][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildDerivativeOperator(d, fes_), fes_);
			for (auto d2{ X }; d2 <= Z; d2++) {
				for (auto f2 : { E, H }) {
					MFN_[f][f2][d] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildOneNormalOperator(f2, { d }, model_, fes_, opts_), fes_);
					MFNN_[f][f2][d][d2] = buildByMult(*buildInverseMassMatrix(f, model_, fes_), *buildTwoNormalOperator(f2, { d, d2 }, model_, fes_, opts_), fes_);
				}
			}
		}
	}

	}
}