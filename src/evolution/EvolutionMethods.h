#pragma once

#include <mfem.hpp>

#include "components/Model.h"
#include "EvolutionOptions.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"
#include "mfemExtension/LinearIntegrators.h"
#include "mfemExtension/LinearForm_IBFI.hpp"

#include "math/EigenMfemTools.h"

namespace maxwell {

	using namespace mfem;
	//using FiniteElementIBFIOperator = std::unique_ptr<mfemExtension::BilinearFormIBFI>;
	using FiniteElementOperator = std::unique_ptr<mfemExtension::BilinearForm>;
	using FiniteElementVector = std::unique_ptr<mfemExtension::LinearFormIBFI>;

	FiniteElementOperator buildByMult(const BilinearForm& op1, const BilinearForm& op2, FiniteElementSpace&);
	//FiniteElementIBFIOperator buildIBFIByMult(const BilinearForm& op1, const mfemExtension::BilinearFormIBFI& op2, FiniteElementSpace& fes);

	FiniteElementOperator buildInverseMassMatrix(const FieldType&, const Model&, FiniteElementSpace&);
	FiniteElementOperator buildDerivativeOperator(const Direction&, FiniteElementSpace&);
	FiniteElementOperator buildFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildCenteredFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildConductivityOperator(const Model& model, FiniteElementSpace& fes);

	FiniteElementOperator buildFluxFunctionOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildFluxOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildPenaltyFixOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildZeroNormalOperator(const FieldType&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildOneNormalOperator (const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildTwoNormalOperator (const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildZeroNormalIBFIOperator(const FieldType& f, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildOneNormalIBFIOperator (const FieldType& f, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildTwoNormalIBFIOperator (const FieldType& f, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildTFSFOperator(const FieldType& f, FiniteElementSpace& fes, double coeff);
	FiniteElementOperator buildOneNormalTotalFieldOperator (const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	std::map<BdrCond, std::vector<double>> bdrCoeffCheck(const EvolutionOptions&);

	FieldType altField(const FieldType& f);

	void calculateEigenvalues(const Eigen::MatrixXd&, Eigen::VectorXcd&);
	void calculateEigenvalues(SparseMatrix& mat, Vector& res);
	void exportSparseToMarketFile(const Eigen::MatrixXd&);
	double usePowerMethod(const Eigen::SparseMatrix<double>&, int iterations);


	std::vector<int> calcOffsetCoeff1D(const std::vector<FieldType>& f);
	std::vector<int> calcOffsetCoeff(const std::vector<FieldType>&, const std::vector<Direction>&);

	void allocateDenseInEigen1D(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const double sign = 1.0);
	void allocateDenseInEigen(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign = 1.0);

}
