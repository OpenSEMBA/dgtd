#pragma once

#include "mfem.hpp"

#include "Types.h"
#include "Model.h"
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
	FiniteElementOperator buildFluxOperator(const FieldType&, const std::vector<Direction>&, bool usePenaltyCoefficients, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
	FiniteElementOperator buildFluxJumpOperator(const FieldType&, const std::vector<Direction>&, bool usePenaltyCoefficients, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
	FiniteElementOperator buildFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
	FiniteElementOperator buildCenteredFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	FiniteElementIBFIOperator buildFluxFunctionOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	FiniteElementOperator buildFluxOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	FiniteElementIBFIOperator buildFluxFunctionOperator1D(Model&, FiniteElementSpace&);
	FiniteElementIBFIOperator buildPenaltyFunctionOperator1D(Model&, FiniteElementSpace&);

	FiniteElementIBFIOperator buildPenaltyIBFIOperator(const FieldType& f, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);
	FiniteElementIBFIOperator buildFluxIBFIOperator(const FieldType& f, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	TFSFOrientationCoefficient interiorBoundaryFaceCoefficient(const BdrCond&);

	FiniteElementOperator buildPenaltyFixOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	FiniteElementOperator buildZeroNormalOperator(const FieldType& f, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);
	FiniteElementOperator buildOneNormalOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);
	FiniteElementOperator buildTwoNormalOperator(const FieldType& f, const std::vector<Direction>& dirTerms, Model& model, FiniteElementSpace& fes, const MaxwellEvolOptions& opts);

	std::map<BdrCond, std::vector<double>> bdrCoeffCheck(const MaxwellEvolOptions&);

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
