#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <../../maxwell/src/maxwell/Types.h>
#include <../../maxwell/src/maxwell/BilinearIntegrators.h>

#include <Eigen/Dense>
#include <maxwell/Model.h>

using namespace mfem;

struct FluxCoefficient {
	double alpha;
	double beta;
};

std::unique_ptr<FiniteElementSpace> buildFiniteElementSpace(const int order)
{
	Mesh mesh = Mesh::MakeCartesian1D(1);
	std::unique_ptr<FiniteElementCollection> fec = std::make_unique<DG_FECollection>(order, 1, BasisType::GaussLobatto);
	std::unique_ptr<FiniteElementSpace> fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get());
	
	return fes;
}

Eigen::MatrixXd convertMFEMDenseToEigen(DenseMatrix* mat)
{
	auto res = Eigen::MatrixXd(mat->Width(), mat->Height());
	for (int i = 0; i < mat->Width(); i++) {
		for (int j = 0; j < mat->Height(); j++) {
			res(i, j) = mat->Elem(i, j);
		}
	}
	return res;
}

Eigen::MatrixXd buildMassMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(fes.get());
	res.AddDomainIntegrator(new MassIntegrator(one));
	res.Assemble();
	res.Finalize();

	return convertMFEMDenseToEigen(res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildInverseMassMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(fes.get());
	res.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	res.Assemble();
	res.Finalize();

	return convertMFEMDenseToEigen(res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildStiffnessMatrixEigen(std::unique_ptr<FiniteElementSpace>& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(fes.get());
	res.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
	res.Assemble();
	res.Finalize();

	return convertMFEMDenseToEigen(res.SpMat().ToDenseMatrix());
}

Eigen::MatrixXd	buildNormalPECFluxOperator1D(std::unique_ptr<FiniteElementSpace>& fes, std::vector<maxwell::Direction> dirVec)
{
	std::vector<maxwell::Direction> dirs = dirVec;
	maxwell::AttributeToBoundary attBdr{ {1,maxwell::BdrCond::PEC},{2,maxwell::BdrCond::PEC} };

	auto res = std::make_unique<BilinearForm>(fes.get());
	{
		FluxCoefficient c = FluxCoefficient{ 0.0, -0.5 };
		res->AddInteriorFaceIntegrator(new maxwell::MaxwellDGTraceJumpIntegrator(dirs, c.beta));
		//res->AddInteriorFaceIntegrator(new DGTraceIntegrator(*(new VectorConstantCoefficient(Vector(1.0))), 0.0, 0.5));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(fes.get()->GetMesh()->bdr_attributes.Max());
	for (auto const& kv : attBdr) {
		Array<int> bdrMarker(fes.get()->GetMesh()->bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int)kv.first - 1] = 1;

		bdrMarkers[(int)kv.first - 1] = bdrMarker;
		FluxCoefficient c = FluxCoefficient{ 0.0, -2.0 };
		res->AddBdrFaceIntegrator(new maxwell::MaxwellDGTraceJumpIntegrator(dirs, c.beta), bdrMarkers[kv.first - 1]);
		//res->AddBdrFaceIntegrator(new DGTraceIntegrator(*(new VectorConstantCoefficient(Vector(1.0))), 0.0, 2.0), bdrMarkers[kv.first - 1]);

		res->Assemble();
		res->Finalize();

		return convertMFEMDenseToEigen(res->SpMat().ToDenseMatrix());
	}
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
	std::unique_ptr<FiniteElementSpace>& fes,
	std::pair<double, double> ab)
{
	auto DGmat = std::make_unique<BilinearForm>(fes.get());
	std::vector<VectorConstantCoefficient> n{ VectorConstantCoefficient(Vector({1.0})) };
	DGmat->AddInteriorFaceIntegrator(
		new DGTraceIntegrator(n[0], ab.first, ab.second));
	DGmat->Assemble();
	DGmat->Finalize();

	std::cout << DGmat.get()->SpMat().ToDenseMatrix() << std::endl;

	return convertMFEMDenseToEigen(DGmat.get()->SpMat().ToDenseMatrix());
}

Eigen::MatrixXd buildEigenMaxwellDGTrace1D(
	std::unique_ptr<FiniteElementSpace>& fes,
	std::vector<maxwell::Direction> dir,
	const double beta)
{
	auto DGmat = std::make_unique<BilinearForm>(fes.get());
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

void checkDenseMatrixSubtractIsValueForAllElem(
	const double val,
	std::unique_ptr<DenseMatrix> m1,
	std::unique_ptr<DenseMatrix> m2)
{
	for (int i = 0; i < m1->Width(); i++) {
		for (int j = 0; j < m1->Height(); j++) {
			EXPECT_NEAR(0.0, m1->Elem(i, j) - m2->Elem(i, j), 1e-3);
		}
	}
}