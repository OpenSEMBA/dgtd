#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "Types.h"
#include "Model.h"
#include "Sources.h"
#include "SourcesManager.h"
#include "EvolutionMethods.h"

namespace maxwell {

class MaxwellEvolution3D_Spectral : public mfem::TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	MaxwellEvolution3D_Spectral(mfem::FiniteElementSpace&, Model&, SourcesManager&, MaxwellEvolOptions&);
	virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

private:

	Eigen::SparseMatrix<double> global_;
	Eigen::SparseMatrix<double> forcing_;
	Eigen::VectorXcd eigenvals_;
	double pmEigenvalue_;

	mfem::FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	MaxwellEvolOptions& opts_;
	

};

}