#pragma once

#include <Eigen/Dense>
#include "evolution/EvolutionMethods.h"

std::unique_ptr<mfem::DenseMatrix> toUnique(mfem::DenseMatrix*);

Eigen::MatrixXd buildMassMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd buildInverseMassMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd build1DStiffnessMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd buildNormalStiffnessMatrixEigen(const maxwell::Direction d, mfem::FiniteElementSpace& fes);

Eigen::MatrixXd	buildNormalSMAFluxOperator1D(
	mfem::FiniteElementSpace&, const std::vector<maxwell::Direction>& dirVec);
Eigen::MatrixXd	buildSMAPenaltyOperator1D(
	mfem::FiniteElementSpace& fes);
