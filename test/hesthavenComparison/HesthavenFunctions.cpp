#include <gtest/gtest.h>

#include "HesthavenFunctions.h"
#include "math/EigenMfemTools.h"
#include "mfemExtension/BilinearIntegrators.h"
#include "components/Model.h"

using namespace mfem;
using namespace maxwell::mfemExtension;

namespace maxwell {

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
	Eigen::Matrix3d blockMat {
		{ 12.0, -16.0,   4.0},
		{  4.0,   0.0,  -4.0},
		{ -4.0,  16.0, -12.0}
	};
	for (int i = 0; i < res.cols(); i += 3) {
		res.block<3, 3>(i, i) = blockMat;
	}
	return res;
}

Eigen::MatrixXd buildMassMatrixEigen(FiniteElementSpace& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(&fes);
	res.AddDomainIntegrator(new MassIntegrator(one));
	res.Assemble();
	res.Finalize();

	return toEigen(*res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildInverseMassMatrixEigen(FiniteElementSpace& fes)
{
	ConstantCoefficient two(2.0);
	BilinearForm res(&fes);
	res.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(two)));
	res.Assemble();
	res.Finalize();

	return toEigen(*res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd build1DStiffnessMatrixEigen(FiniteElementSpace& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(&fes);
	res.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
	res.Assemble();
	res.Finalize();

	return toEigen(*res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildNormalStiffnessMatrixEigen(const Direction d, FiniteElementSpace& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(&fes);
	res.AddDomainIntegrator(new DerivativeIntegrator(one, d));
	res.Assemble();
	res.Finalize();

	return toEigen(*res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd	buildNormalSMAFluxOperator1D(
	FiniteElementSpace& fes, const std::vector<Direction>& dirVec)
{
	std::vector<Direction> dirs = dirVec;
	GeomTagToBoundary attBdr{ {1,BdrCond::SMA},{2,BdrCond::SMA} };
	VectorConstantCoefficient one(Vector({ 1.0 }));

	BilinearForm res(&fes);
	{
		res.AddInteriorFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirVec, 1.0));
		//res.AddInteriorFaceIntegrator(new DGTraceIntegrator(one, 0.0, -1.0));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(fes.GetMesh()->bdr_attributes.Max());
	for (auto const& kv : attBdr) {
		Array<int> bdrMarker(fes.GetMesh()->bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int)kv.first - 1] = 1;
		bdrMarkers[(int)kv.first - 1] = bdrMarker;
		res.AddBdrFaceIntegrator(
			new mfemExtension::MaxwellDGTraceJumpIntegrator(dirVec, 1.0), bdrMarkers[kv.first - 1]);
		//res.AddBdrFaceIntegrator(new DGTraceIntegrator(one, 0.0, -1.0),bdrMarkers[kv.first -1]);

	}
	res.Assemble();
	res.Finalize();

	return toEigen(*res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd	buildSMAPenaltyOperator1D(
	FiniteElementSpace& fes)
{
	GeomTagToBoundary attBdr{ {1,BdrCond::SMA},{2,BdrCond::SMA} };
	VectorConstantCoefficient one(Vector({ 1.0 }));

	BilinearForm res(&fes);
	{
		res.AddInteriorFaceIntegrator(new DGTraceIntegrator(one, 0.0, 1.0));
		//res.AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator(std::vector<Direction>{}, -1.0));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(fes.GetMesh()->bdr_attributes.Max());
	for (auto const& kv : attBdr) {
		Array<int> bdrMarker(fes.GetMesh()->bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int)kv.first - 1] = 1;
		bdrMarkers[(int)kv.first - 1] = bdrMarker;
		res.AddBdrFaceIntegrator(new DGTraceIntegrator(one, 0.0, 1.0), bdrMarkers[kv.first - 1]);
		//res.AddBdrFaceIntegrator(new MaxwellDGTraceJumpIntegrator(std::vector<Direction>{}, -1.0), bdrMarkers[kv.first - 1]);

	}
	res.Assemble();
	res.Finalize();

	return toEigen(*res.SpMat().ToDenseMatrix());
}

}