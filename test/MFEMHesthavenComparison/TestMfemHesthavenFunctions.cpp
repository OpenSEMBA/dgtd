#include "TestMfemHesthavenFunctions.h"
#include "../TestGlobalFunctions.h"
#include "mfem.hpp"
#include <Eigen/Dense>

using namespace mfem;


std::unique_ptr<FiniteElementSpace> buildFiniteElementSpace(const int order)
{
	Mesh mesh = Mesh::MakeCartesian1D(1);
	std::unique_ptr<FiniteElementCollection> fec = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get());
	
	return fes;
}


Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(
	const int order,
	const int elements
)
{
	std::unique_ptr<DenseMatrix> res = std::make_unique<DenseMatrix>((order + 1) * elements);
	res->operator=(0.0);

	for (int i = 1; i <= order; i++) {
		for (int j = 1; j <= order; j++) {
			for (int it = 0; it < elements - 1; it++) {
				if (i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == order + ((order + 1) * it) ||
					i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == (order + 1) + ((order + 1) * it)) {
					res->Elem(i + (order + 1) * it, j + (order + 1) * it) = 0.5;
					res->Elem(i + (order + 1) * it, j + 1 + (order + 1) * it) = 0.5;
					res->Elem(i + 1 + (order + 1) * it, j + (order + 1) * it) = -0.5;
					res->Elem(i + 1 + (order + 1) * it, j + 1 + (order + 1) * it) = -0.5;
				}
			}
		}
	}
	return convertMFEMDenseToEigen(res.get());
}

Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(
	const int order,
	const int elements
)
{
	std::unique_ptr<DenseMatrix> res = std::make_unique<DenseMatrix>((order + 1) * elements);
	res->operator=(0.0);

	for (int i = 1; i <= order; i++) {
		for (int j = 1; j <= order; j++) {
			for (int it = 0; it < elements - 1; it++) {
				if (i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == order + ((order + 1) * it) ||
					i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == (order + 1) + ((order + 1) * it)) {
					res->Elem(i + (order + 1) * it, j + (order + 1) * it) = 1.0;
					res->Elem(i + (order + 1) * it, j + 1 + (order + 1) * it) = -1.0;
					res->Elem(i + 1 + (order + 1) * it, j + (order + 1) * it) = -1.0;
					res->Elem(i + 1 + (order + 1) * it, j + 1 + (order + 1) * it) = 1.0;
				}
			}
		}
	}
	return convertMFEMDenseToEigen(res.get());
}

Eigen::MatrixXd buildEigenDGTrace1D(
	FiniteElementSpace* fes,
	std::pair<double, double> ab)
{
	auto DGmat = std::make_unique<BilinearForm>(fes);
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat->AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], ab.first, ab.second));
	DGmat->Assemble();
	DGmat->Finalize();

	std::cout << DGmat.get()->SpMat().ToDenseMatrix() << std::endl;

	return convertMFEMDenseToEigen(DGmat.get()->SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildEigenMaxwellDGTrace1D(
	FiniteElementSpace* fes,
	std::vector<maxwell::Direction> dir,
	const double beta)
{
	auto DGmat = std::make_unique<BilinearForm>(fes);
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat->AddInteriorFaceIntegrator(
		new maxwell::MaxwellDGTraceJumpIntegrator(dir, beta));
	DGmat->Assemble();
	DGmat->Finalize();

	std::cout << DGmat.get()->SpMat().ToDenseMatrix() << std::endl;

	return convertMFEMDenseToEigen(DGmat.get()->SpMat().ToDenseMatrix());
}

std::unique_ptr<BilinearForm> buildBilinearFormWith1DCartesianMesh(
	const int elements,
	const int order,
	std::pair<double, double> ab)
{
	Mesh mesh = Mesh::MakeCartesian1D(elements);
	FiniteElementCollection* fec = new DG_FECollection(order, mesh.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);
	auto DGmat = std::make_unique<BilinearForm>(fes);
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat->AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], ab.first, ab.second));
	DGmat->Assemble();
	DGmat->Finalize();
	return DGmat;
}

std::unique_ptr<BilinearForm> buildMaxwellBilinearFormWith1DCartesianMesh(
	const int elements,
	const int order,
	std::vector<maxwell::Direction> dir,
	const double beta)
{
	Mesh mesh = Mesh::MakeCartesian1D(elements);
	FiniteElementCollection* fec = new DG_FECollection(order, mesh.Dimension(), BasisType::GaussLobatto);
	FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);
	auto DGmat = std::make_unique<BilinearForm>(fes);
	DGmat->AddInteriorFaceIntegrator(
		new maxwell::MaxwellDGTraceJumpIntegrator(dir, beta));
	DGmat->Assemble();
	DGmat->Finalize();
	return DGmat;

}

Eigen::Matrix<double, 27, 27> build3DOneElementDMatrix()
{
	auto res = Eigen::Matrix<double, 27, 27>();
	res.setZero();
	auto blockMat = Eigen::Matrix3d{
			{-1.5, 2.0,-0.5},
			{-0.5, 0.0, 0.5},
			{ 0.5,-2.0, 1.5} };
	for (int i = 0; i < res.cols(); i += 3) {
		res.block<3, 3>(i, i) = blockMat;
	}
	return res;
}