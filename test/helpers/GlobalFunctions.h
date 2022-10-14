#pragma once

#include <Eigen/Dense>

#include <mfem.hpp>
#include "maxwell/Types.h"

std::unique_ptr<mfem::DenseMatrix> toUnique(mfem::DenseMatrix*);
Eigen::MatrixXd toEigen(const mfem::DenseMatrix&);

Eigen::MatrixXd buildMassMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd buildInverseMassMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd buildStiffnessMatrixEigen(mfem::FiniteElementSpace&);

Eigen::MatrixXd	buildNormalSMAFluxOperator1D(
	mfem::FiniteElementSpace&, const std::vector<maxwell::Direction>& dirVec);
Eigen::MatrixXd	buildSMAPenaltyOperator1D(
	mfem::FiniteElementSpace& fes);