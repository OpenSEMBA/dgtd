#include <Eigen/Dense>
#include "mfem.hpp"
#include "maxwell/Types.h"
#include "maxwell/Model.h"
#include "maxwell/BilinearIntegrators.h"

using namespace mfem;

struct FluxCoefficient {
	double alpha;
	double beta;
};

Eigen::MatrixXd convertMFEMDenseToEigen(const DenseMatrix* mat);
Eigen::MatrixXd buildMassMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes);
Eigen::MatrixXd buildInverseMassMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes);
Eigen::MatrixXd buildStiffnessMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes);
Eigen::MatrixXd	buildNormalPECFluxOperator1D(FiniteElementSpace* fes, std::vector<maxwell::Direction> dirVec);