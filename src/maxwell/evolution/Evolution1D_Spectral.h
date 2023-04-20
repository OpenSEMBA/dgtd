#pragma once

#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/LinearIntegrators.h"

#include "Types.h"
#include "Model.h"
#include "Sources.h"
#include "SourcesManager.h"
#include "EvolutionMethods.h"

namespace maxwell {

class MaxwellEvolution1D_Spectral : public mfem::TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 1;

	MaxwellEvolution1D_Spectral(mfem::FiniteElementSpace&, Model&, SourcesManager&, MaxwellEvolOptions&);
	virtual void Mult(const Vector& x, Vector& y) const;
	double GetTime() const { return t; }
	void SetTime(const double time) { t = time; }

	const mfem::FiniteElementSpace& getFES() { return fes_; }

private:

	Eigen::SparseMatrix<double> global_;
	Eigen::SparseMatrix<double> forcing_;
	Vector eigenvals_;
	double pmEigenvalue_;

	mfem::FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	MaxwellEvolOptions& opts_;

};

}