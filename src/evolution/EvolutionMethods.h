#pragma once

#include <mfem.hpp>

#include "components/Model.h"
#include "EvolutionOptions.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "mfemExtension/BilinearForm_IBFI.hpp"
#include "mfemExtension/LinearIntegrators.h"
#include "mfemExtension/LinearForm_IBFI.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

namespace maxwell {

	using namespace mfem;
	using FiniteElementIBFIOperator = std::unique_ptr<mfemExtension::BilinearFormIBFI>;
	using FiniteElementOperator = std::unique_ptr<mfemExtension::BilinearForm>;
	using FiniteElementVector = std::unique_ptr<mfemExtension::LinearFormIBFI>;


	Eigen::MatrixXd toEigen(const DenseMatrix&);

	FiniteElementOperator buildByMult(const BilinearForm& op1, const BilinearForm& op2, FiniteElementSpace&);
	FiniteElementIBFIOperator buildIBFIByMult(const BilinearForm& op1, const mfemExtension::BilinearFormIBFI& op2, FiniteElementSpace& fes);

	FiniteElementOperator buildInverseMassMatrix(const FieldType&, const Model&, FiniteElementSpace&);
	FiniteElementOperator buildDerivativeOperator(const Direction&, FiniteElementSpace&);
	FiniteElementOperator buildFluxOperator(const FieldType&, const std::vector<Direction>&, bool usePenaltyCoefficients, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildFluxJumpOperator(const FieldType&, const std::vector<Direction>&, bool usePenaltyCoefficients, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildCenteredFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementIBFIOperator buildFluxFunctionOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildFluxOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementIBFIOperator buildFluxFunctionOperator1D(Model&, FiniteElementSpace&);
	FiniteElementIBFIOperator buildPenaltyFunctionOperator1D(Model&, FiniteElementSpace&);

	TFSFOrientationCoefficient interiorBoundaryFaceCoefficient(const BdrCond&);

	FiniteElementOperator buildPenaltyFixOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);

	FiniteElementOperator buildZeroNormalOperator(const FieldType&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildOneNormalOperator (const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementOperator buildTwoNormalOperator (const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementIBFIOperator buildZeroNormalIBFIOperator(const FieldType& f, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementIBFIOperator buildOneNormalIBFIOperator (const FieldType& f, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);
	FiniteElementIBFIOperator buildTwoNormalIBFIOperator (const FieldType& f, const std::vector<Direction>&, Model&, FiniteElementSpace&, const EvolutionOptions&);


	std::map<BdrCond, std::vector<double>> bdrCoeffCheck(const EvolutionOptions&);

	FieldType altField(const FieldType& f);

	void calculateEigenvalues(const Eigen::MatrixXd&, Eigen::VectorXcd&);
	void calculateEigenvalues(SparseMatrix& mat, Vector& res);
	void exportSparseToMarketFile(const Eigen::MatrixXd&);
	double usePowerMethod(const Eigen::SparseMatrix<double>&, int iterations);

	Eigen::VectorXd toEigenVector(const Vector&);
	Eigen::VectorXcd toComplexEigenVector(const Vector&);
	Vector toMFEMVector(const Eigen::VectorXd&);
	SparseMatrix toMFEMSparse(const Eigen::SparseMatrix<double>&);

	std::vector<int> calcOffsetCoeff1D(const std::vector<FieldType>& f);
	std::vector<int> calcOffsetCoeff(const std::vector<FieldType>&, const std::vector<Direction>&);


	void allocateDenseInEigen1D(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const double sign = 1.0);
	void allocateDenseInEigen(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign = 1.0);

}