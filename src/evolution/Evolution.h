#pragma once

#include "mfemExtension/BilinearIntegrators.h"

#include "components/Model.h"
#include "solver/SourcesManager.h"
#include "EvolutionMethods.h"
#include "components/SubMesher.h"

#include "evolution/EvolutionMethods.h"
#include "evolution/HesthavenEvolutionTools.h"

#include <chrono>

namespace maxwell {

using MatrixStorageLT = std::set<DynamicMatrix, MatrixCompareLessThan>;

class Evolution: public mfem::TimeDependentOperator {
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	Evolution(mfem::FiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
	virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

private:

	std::array<FiniteElementOperator, 2> MInv_;
	std::array<FiniteElementOperator, 2> MInvTFSF_;

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_;
	std::array<FiniteElementOperator, 2> MP_;
	std::array<FiniteElementOperator, 2> MFF_;

	Vector CND_;

	std::array<std::array<FiniteElementOperator, 3>,2> MBF_;

	std::array<FiniteElementOperator, 2> MPB_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFNB_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNNB_;

	std::array<FiniteElementOperator, 2> MTF_;
	std::array<FiniteElementOperator, 2> MSF_;

	//Total Field and Scattered Field operators for SubMeshing

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_TF_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_TF_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_TF_;
	std::array<FiniteElementOperator, 2> MP_TF_;

	std::array<std::array<FiniteElementOperator, 3>, 2> MS_SF_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_SF_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_SF_;
	std::array<FiniteElementOperator, 2> MP_SF_;


	std::array<std::array<FiniteElementOperator, 3>, 2> MS_GTFSF_;
	std::array<std::array<std::array<std::array<FiniteElementOperator, 3>, 3>, 2>, 2> MFNN_GTFSF_;
	std::array<std::array<std::array<FiniteElementOperator, 3>, 2>, 2> MFN_GTFSF_;
	std::array<FiniteElementOperator, 2> MP_GTFSF_;

	/* */

	mfem::FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;
	

};

class HesthavenEvolution : public mfem::TimeDependentOperator 
{
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	HesthavenEvolution(mfem::FiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
	virtual void Mult(const mfem::Vector& in, mfem::Vector& out) const;

private:

	void assembleTFSFConnectivity(const DynamicMatrix& matrix, FaceElementTransformations*, double faceOri);
	void evaluateTFSF(HesthavenFields& jumps) const;
	void initBdrConnectivityMaps(const std::vector<std::vector<NodeId>>& bdr2nodes);
	const Eigen::VectorXd applyScalingFactors(const ElementId, const Eigen::VectorXd& flux) const;
	FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;
	MatrixStorageLT matrixStorage_;
	std::vector<HesthavenElement> hestElemStorage_;
	DynamicMatrix refInvMass_;
	DynamicMatrix refLIFT_;
	GlobalConnectivity connectivity_;
	BoundaryMaps bdr_connectivity_;
	InteriorBoundaryMaps int_bdr_connectivity_;
	TFSFConnectivity tfsf_connectivity_;

};

class CurvedEvolution : public mfem::TimeDependentOperator 
{
public:
	static const int numberOfFieldComponents = 2;
	static const int numberOfMaxDimensions = 3;

	CurvedEvolution(mfem::FiniteElementSpace&, Model&, SourcesManager&, EvolutionOptions&);
	virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

private:

	mfem::FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;
};

}