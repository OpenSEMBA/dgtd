#pragma once

#include "mfem.hpp"
#include "maxwell/Types.h"
#include "maxwell/BilinearIntegrators.h"

#include <Eigen/Dense>
#include <maxwell/Model.h>

Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(const int order,const int elements);
Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(const int order, const int elements);
Eigen::MatrixXd buildDGTrace1DEigen(FiniteElementSpace& fes, maxwell::FluxCoefficient ab);
Eigen::MatrixXd buildMaxwellDGTrace1DEigen(FiniteElementSpace& fes, const std::vector<maxwell::Direction>& dir, const double beta);
Eigen::MatrixXd build3DOneElementDMatrix();