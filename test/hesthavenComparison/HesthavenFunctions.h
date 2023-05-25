#pragma once

#include <mfem.hpp>
#include <Eigen/Dense>

#include "components/Types.h"

namespace maxwell {
Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(const int order,const int elements);

Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(const int order, const int elements);

Eigen::MatrixXd buildMaxwellDGTrace1DEigen(
	mfem::FiniteElementSpace& fes, const std::vector<maxwell::Direction>& dir, const double beta);

Eigen::MatrixXd build3DOneElementDMatrix();

Eigen::MatrixXd buildMatrixForMSTest4E();

Eigen::MatrixXd buildMassMatrixEigen(mfem::FiniteElementSpace&);

Eigen::MatrixXd buildInverseMassMatrixEigen(mfem::FiniteElementSpace&);

Eigen::MatrixXd build1DStiffnessMatrixEigen(mfem::FiniteElementSpace&);

Eigen::MatrixXd buildNormalStiffnessMatrixEigen(
	const maxwell::Direction d, mfem::FiniteElementSpace& fes);

Eigen::MatrixXd	buildNormalSMAFluxOperator1D(
	mfem::FiniteElementSpace&, const std::vector<maxwell::Direction>& dirVec);

Eigen::MatrixXd	buildSMAPenaltyOperator1D(
	mfem::FiniteElementSpace& fes);

}