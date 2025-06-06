#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "components/Model.h"
#include "components/SubMesher.h"

#include "evolution/MaxwellEvolutionMethods.h"
#include "evolution/HesthavenEvolutionMethods.h"

namespace maxwell {

	using straightElemList = mfem::Array<ElementId>;
	using curvedElemMap = std::map<ElementId, mfem::Array<NodeId>>;
	std::pair<straightElemList, curvedElemMap> initCurvedAndLinearElementsLists(const mfem::ParFiniteElementSpace& fes, const std::vector<Source::Position>& curved_pos);

class HesthavenEvolution : public mfem::TimeDependentOperator
{
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	HesthavenEvolution(mfem::ParFiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
	virtual void Mult(const mfem::Vector& in, mfem::Vector& out) const;

	const HesthavenElement& getHesthavenElement(const ElementId& e) const { return hestElemLinearStorage_[e]; }

private:

	void evaluateTFSF(HesthavenFields& jumps) const;
	void storeDirectionalMatrices(ParFiniteElementSpace& subFES, const DynamicMatrix& refInvMass, HesthavenElement&);
	void checkForTFSFInCurvedElements();
	void applyBoundaryConditionsToNodes(const BoundaryMaps&, const FieldsInputMaps& in, HesthavenFields& out) const;
	bool isDoFinCurvedElement(const NodeId& d) const;
	
	const Eigen::VectorXd applyLIFT(const Eigen::VectorXd& fscale, Eigen::VectorXd& flux) const;
	
	ParFiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;

	Array<ElementId> linearElements_;
	std::map<ElementId, Array<NodeId>> curvedElements_;

	std::set<DynamicMatrix, MatrixCompareLessThan> matrixStorage_;
	std::vector<HesthavenElement> hestElemLinearStorage_;
	std::vector<HesthavenCurvedElement> hestElemCurvedStorage_;
	DynamicMatrix refLIFT_;
	Connectivities connectivity_;
	std::vector<Source::Position> positions_;

};
}