#include <gtest/gtest.h>

#include "maxwell/mfemExtension/BilinearIntegrators.h"

#include "MfemHesthavenFunctionsTest.h"
#include "GlobalFunctions.h"

using namespace mfem;
using namespace maxwell::mfemExtension;

Eigen::MatrixXd buildExpectedAverageDenseMatrix1D(
	const int order,
	const int elements
)
{
	DenseMatrix res((order + 1) * elements);
	res = 0.0;

	for (int i = 1; i <= order; i++) {
		for (int j = 1; j <= order; j++) {
			for (int it = 0; it < elements - 1; it++) {
				if (i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == order + ((order + 1) * it) ||
					i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == (order + 1) + ((order + 1) * it)) {
					res.Elem(i + (order + 1) * it, j + (order + 1) * it) = 0.5;
					res.Elem(i + (order + 1) * it, j + 1 + (order + 1) * it) = 0.5;
					res.Elem(i + 1 + (order + 1) * it, j + (order + 1) * it) = -0.5;
					res.Elem(i + 1 + (order + 1) * it, j + 1 + (order + 1) * it) = -0.5;
				}
			}
		}
	}
	return maxwell::toEigen(res);
}

Eigen::MatrixXd buildExpectedJumpDenseMatrix1D(
	const int order,
	const int elements
)
{
	DenseMatrix res((order + 1) * elements);
	res = 0.0;

	for (int i = 1; i <= order; i++) {
		for (int j = 1; j <= order; j++) {
			for (int it = 0; it < elements - 1; it++) {
				if (i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == order + ((order + 1) * it) ||
					i + (order + 1) * it == order + ((order + 1) * it) && j + (order + 1) * it == (order + 1) + ((order + 1) * it)) {
					res.Elem(i + (order + 1) * it, j + (order + 1) * it) = 1.0;
					res.Elem(i + (order + 1) * it, j + 1 + (order + 1) * it) = -1.0;
					res.Elem(i + 1 + (order + 1) * it, j + (order + 1) * it) = -1.0;
					res.Elem(i + 1 + (order + 1) * it, j + 1 + (order + 1) * it) = 1.0;
				}
			}
		}
	}
	return maxwell::toEigen(res);
}

Eigen::MatrixXd buildDGTraceAverage1DEigen(
	FiniteElementSpace& fes,
	maxwell::FluxCoefficient ab)
{
	BilinearForm DGmat(&fes);
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], ab.beta, 0.0));
	DGmat.Assemble();
	DGmat.Finalize();

	return maxwell::toEigen(*DGmat.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildDGTraceJump1DEigen(
	FiniteElementSpace& fes,
	maxwell::FluxCoefficient ab)
{
	BilinearForm DGmat(&fes);
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat.AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], 0.0, ab.beta));
	DGmat.Assemble();
	DGmat.Finalize();

	return maxwell::toEigen(*DGmat.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildMaxwellDGTrace1DEigen(
	FiniteElementSpace& fes,
	const std::vector<maxwell::Direction>& dir,
	const double beta)
{
	BilinearForm DGmat(&fes);
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat.AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dir, beta));
	DGmat.Assemble();
	DGmat.Finalize();

	return maxwell::toEigen(*DGmat.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd build3DOneElementDMatrix()
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

Eigen::MatrixXd buildMatrixForMSTest4E()
{
	auto res = Eigen::Matrix<double, 12, 12>();
	res.setZero();
	Eigen::Matrix3d blockMat{
		{ 12.0, -16.0,   4.0},
		{  4.0,   0.0,  -4.0},
		{ -4.0,  16.0, -12.0}};
	for (int i = 0; i < res.cols(); i += 3) {
		res.block<3, 3>(i, i) = blockMat;
	}
	return res;
}