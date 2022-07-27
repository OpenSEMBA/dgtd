#include <Eigen/Dense>
#include "mfem.hpp"

using namespace mfem;

Eigen::MatrixXd convertMFEMDenseToEigen(DenseMatrix* mat);
Eigen::MatrixXd buildMassMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes);
Eigen::MatrixXd buildInverseMassMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes);
Eigen::MatrixXd buildStiffnessMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes);