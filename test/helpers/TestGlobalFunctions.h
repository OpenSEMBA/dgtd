#include <Eigen/Dense>

#include "mfem.hpp"
#include "maxwell/Types.h"

std::unique_ptr<mfem::DenseMatrix> toUnique(mfem::DenseMatrix*);
Eigen::MatrixXd convertMFEMDenseToEigen(const mfem::DenseMatrix&);

Eigen::MatrixXd buildMassMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd buildInverseMassMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd buildStiffnessMatrixEigen(mfem::FiniteElementSpace&);
Eigen::MatrixXd	buildNormalPECFluxOperator1D(mfem::FiniteElementSpace&, std::vector<maxwell::Direction> dirVec);