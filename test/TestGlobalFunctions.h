#include <Eigen/Dense>
#include "mfem.hpp"
#include "maxwell/Types.h"
#include "maxwell/Model.h"
#include "maxwell/BilinearIntegrators.h"

using namespace mfem;

Eigen::MatrixXd convertMFEMDenseToEigen(const DenseMatrix* mat);
Eigen::MatrixXd buildMassMatrixEigen(FiniteElementSpace* fes);
Eigen::MatrixXd buildInverseMassMatrixEigen(FiniteElementSpace* fes);
Eigen::MatrixXd buildStiffnessMatrixEigen(FiniteElementSpace* fes);
Eigen::MatrixXd	buildNormalPECFluxOperator1D(FiniteElementSpace* fes, std::vector<maxwell::Direction> dirVec);