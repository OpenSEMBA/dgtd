#ifndef TESTMFEMHESTHAVENFUNCTIONS_H
#define TESTMFEMHESTHAVENFUNCTIONS_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <../../maxwell/src/maxwell/Types.h>
#include <../../maxwell/src/maxwell/BilinearIntegrators.h>
#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <maxwell/Model.h>


std::unique_ptr<FiniteElementSpace> buildFiniteElementSpace(const int order);
Eigen::MatrixXd	buildNormalPECFluxOperator1D(std::unique_ptr<FiniteElementSpace>& fes, std::vector<maxwell::Direction> dirVec);
Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(const int order,const int elements);
Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(const int order, const int elements);
Eigen::MatrixXd buildEigenDGTrace1D(std::unique_ptr<FiniteElementSpace>& fes, std::pair<double, double> ab);
Eigen::MatrixXd buildEigenMaxwellDGTrace1D(std::unique_ptr<FiniteElementSpace>& fes, std::vector<maxwell::Direction> dir, const double beta);
std::unique_ptr<BilinearForm> buildBilinearFormWith1DCartesianMesh(const int elements, const int order, std::pair<double, double> ab);
std::unique_ptr<BilinearForm> buildMaxwellBilinearFormWith1DCartesianMesh(const int elements, const int order, std::vector<maxwell::Direction> dir, const double beta);
Eigen::Matrix<double, 27, 27> build3DOneElementDMatrix();

#endif