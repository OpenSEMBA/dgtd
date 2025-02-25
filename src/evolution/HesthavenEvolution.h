#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "components/Model.h"
#include "components/SubMesher.h"

#include "evolution/MaxwellEvolutionMethods.h"
#include "evolution/HesthavenEvolutionMethods.h"

namespace maxwell {

class HesthavenEvolution : public mfem::TimeDependentOperator
{
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	HesthavenEvolution(mfem::FiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
	virtual void Mult(const mfem::Vector& in, mfem::Vector& out) const;

	const HesthavenElement& getHesthavenElement(const ElementId& e) const { return hestElemStorage_[e]; }

private:

	void evaluateTFSF(HesthavenFields& jumps) const;
	void storeDirectionalMatrices(FiniteElementSpace& subFES, const DynamicMatrix& refInvMass, HesthavenElement&);
	
	const Eigen::VectorXd applyLIFT(const ElementId, Eigen::VectorXd& flux) const;
	
	FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;

	std::set<DynamicMatrix, MatrixCompareLessThan> matrixStorage_;
	std::vector<HesthavenElement> hestElemStorage_;
	DynamicMatrix refLIFT_;
	Connectivities connectivity_;
	std::vector<Source::Position> positions_;

};
}