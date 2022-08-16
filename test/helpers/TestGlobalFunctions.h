#include <Eigen/Dense>
#include "mfem.hpp"
#include "maxwell/Types.h"
#include "maxwell/Model.h"
#include "maxwell/BilinearIntegrators.h"

using namespace mfem;

std::unique_ptr<DenseMatrix> toUnique(DenseMatrix*);
Eigen::MatrixXd convertMFEMDenseToEigen(const DenseMatrix&);

Eigen::MatrixXd buildMassMatrixEigen(FiniteElementSpace&);
Eigen::MatrixXd buildInverseMassMatrixEigen(FiniteElementSpace&);
Eigen::MatrixXd buildStiffnessMatrixEigen(FiniteElementSpace&);
Eigen::MatrixXd	buildNormalPECFluxOperator1D(FiniteElementSpace&, std::vector<maxwell::Direction> dirVec);