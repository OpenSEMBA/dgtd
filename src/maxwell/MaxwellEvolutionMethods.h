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
	void checkEigenvalues(const Eigen::VectorXcd&);
	void exportSparseToMarketFile(const Eigen::MatrixXd&);

	Eigen::VectorXd toEigenVector(const Vector&);
	Vector toMFEMVector(const Eigen::VectorXd&);

	std::vector<int> calcOffsetCoeff(const std::vector<FieldType>&, const std::vector<Direction>&);

template <class T>
void allocateDenseInEigen1D(const std::array<T, 2>& arr, Eigen::MatrixXd& res, const double sign = 1.0, bool altField = false)
{
	int offset = arr[E]->SpMat().ToDenseMatrix()->Height();
	for (int i = 0; i < arr[E]->Height(); ++i) {
		for (int j = 0; j < arr[E]->Width(); ++j) {
			switch (altField) {
			case false:
				res(i, j) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j + offset) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			case true:
				res(i, j + offset) += sign * arr[E]->SpMat().ToDenseMatrix()->Elem(i, j);
				res(i + offset, j) += sign * arr[H]->SpMat().ToDenseMatrix()->Elem(i, j);
				break;
			}
		}
	}
}


template <class T>
void allocateDenseInEigen(const std::unique_ptr<T>& bilForm, Eigen::MatrixXd& res, const std::vector<FieldType> f, const std::vector<Direction> d, const double sign = 1.0)
{
	FiniteElementOperator F;
	bilForm->SpMat().ToDenseMatrix()->Print(std::cout);
	std::cout << "AAAAAAAAAAAAAAA" << std::endl;
	int offset = bilForm->SpMat().ToDenseMatrix()->Height();
	std::vector<int> offsetCoeff{ calcOffsetCoeff(f,d) };

	for (int i = 0; i < bilForm->SpMat().ToDenseMatrix()->Height(); ++i) {
		for (int j = 0; j < bilForm->SpMat().ToDenseMatrix()->Width(); ++j) {
			if (bilForm->SpMat().ToDenseMatrix()->Elem(i, j) != 0.0) {
				res(i + offset * offsetCoeff[0], j + offset * offsetCoeff[1]) += sign * bilForm->SpMat().ToDenseMatrix()->Elem(i, j);
			}
		}
	}
}

}
