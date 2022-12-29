#pragma once

#include <mfem.hpp>
#include "maxwell/Types.h"

#include <Eigen/Dense>
#include <maxwell/Model.h>

Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(const int order,const int elements);
Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(const int order, const int elements);

Eigen::MatrixXd buildDGTraceAverage1DEigen(
	mfem::FiniteElementSpace& fes, maxwell::FluxCoefficient ab);
Eigen::MatrixXd buildDGTraceJump1DEigen(
	mfem::FiniteElementSpace& fes, maxwell::FluxCoefficient ab);
Eigen::MatrixXd buildMaxwellDGTrace1DEigen(
	mfem::FiniteElementSpace& fes, const std::vector<maxwell::Direction>& dir, const double beta);

Eigen::MatrixXd build3DOneElementDMatrix();

Eigen::MatrixXd buildMatrixForMSTest4E();
