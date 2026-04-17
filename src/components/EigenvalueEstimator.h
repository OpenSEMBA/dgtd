#pragma once

#include "evolution/MaxwellEvolutionMethods.h"
#include "solver/Solver.h"

namespace maxwell {

class EigenvalueEstimator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	EigenvalueEstimator(mfem::ParFiniteElementSpace&, Model&, EvolutionOptions&);

	const Eigen::MatrixXd& getElementMatrix() { return mat_; };

private:

	int getOffset(const FieldType&, const Direction&);

	mfem::ParFiniteElementSpace& fes_;
	Model& model_;
	EvolutionOptions& opts_;

	Eigen::MatrixXd mat_;

	std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 2> MS_;
	std::array<std::array<std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 3>, 2>, 2> MFNN_;
	std::array<std::array<std::array<std::unique_ptr<ParBilinearForm>, 3>, 2>, 2> MFN_;
	std::array<std::unique_ptr<ParBilinearForm>, 2> MP_;

};
}