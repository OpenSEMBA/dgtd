#pragma once

#include <Eigen/Dense>
#include "components/Types.h"

Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(const int order,const int elements);
Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(const int order, const int elements);

Eigen::MatrixXd buildMaxwellDGTrace1DEigen(
	mfem::FiniteElementSpace& fes, const std::vector<maxwell::Direction>& dir, const double beta);

Eigen::MatrixXd build3DOneElementDMatrix();

Eigen::MatrixXd buildMatrixForMSTest4E();
