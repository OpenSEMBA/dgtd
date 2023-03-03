#pragma once

#include "Types.h"
#include "mfem.hpp"
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
	FiniteElementOperator buildFluxOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	FiniteElementIBFIOperator buildFluxFunctionOperator(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementIBFIOperator buildPenaltyFunctionOperator(const FieldType&, Model&, FiniteElementSpace&);

	FiniteElementOperator buildFluxOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&);
	FiniteElementOperator buildPenaltyOperator1D(const FieldType&, const std::vector<Direction>&, Model&, FiniteElementSpace&, const MaxwellEvolOptions&);

	FiniteElementIBFIOperator buildFluxFunctionOperator1D(Model&, FiniteElementSpace&);
	FiniteElementIBFIOperator buildPenaltyFunctionOperator1D(Model&, FiniteElementSpace&);

	TFSFOrientationCoefficient interiorBoundaryFaceCoefficient(const BdrCond&);

	FluxCoefficient interiorFluxCoefficient();
	FluxCoefficient interiorPenaltyFluxCoefficient(const MaxwellEvolOptions&);
	FluxCoefficient boundaryFluxCoefficient(const FieldType&, const BdrCond&);
	FluxCoefficient boundaryPenaltyFluxCoefficient(const FieldType&, const BdrCond&, const MaxwellEvolOptions&);

	FieldType altField(const FieldType& f);

	void calculateEigenvalues(const Eigen::MatrixXd&, Eigen::VectorXcd&);
	void calculateEigenvalues(SparseMatrix& mat, Vector& res);
	void checkEigenvalues(const Eigen::VectorXcd&);
	void exportSparseToMarketFile(const Eigen::MatrixXd&);
	double usePowerMethod(const Eigen::SparseMatrix<double>&, int iterations);

	Eigen::VectorXd toEigenVector(const Vector&);
	Vector toMFEMVector(const Eigen::VectorXd&);
	SparseMatrix toMFEMSparse(const Eigen::SparseMatrix<double>&);

	std::vector<int> calcOffsetCoeff1D(const std::vector<FieldType>& f);
	std::vector<int> calcOffsetCoeff(const std::vector<FieldType>&, const std::vector<Direction>&);

	void allocateDenseInEigen1D(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const double sign = 1.0);
	void allocateDenseInEigen(DenseMatrix* bilForm, Eigen::SparseMatrix<double>& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign = 1.0);

}
