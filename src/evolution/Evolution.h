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
	void Mult(const mfem::Vector& in, mfem::Vector& out);

private:

	mfem::FiniteElementSpace& fes_;
	Model& model_;
	SourcesManager& srcmngr_;
	EvolutionOptions& opts_;
	std::set<DynamicMatrix, MatrixCompareLessThan> matrixStorage_;
	std::vector<HesthavenElement> hestElemStorage_;
	GlobalConnectivityMap connectivity_;
	GlobalBoundaryMap bdr_connectivity_;
	GlobalInteriorBoundaryMap int_bdr_connectivity_;

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