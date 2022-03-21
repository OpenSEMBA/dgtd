#include "gtest/gtest.h"

#include "maxwell1D/Solver.h"
#include "maxwell1D/BilinearIntegrators.h"

using namespace Maxwell1D;

class TestMaxwellDGTrace : public ::testing::Test {

};

TEST_F(TestMaxwellDGTrace, DGTraces)
{
	const int order = 1;
	const int dimension = 1;
	const double tol = 1e-3;
	Mesh mesh = Mesh::MakeCartesian1D(10);
	FiniteElementCollection* fec = new DG_FECollection(order, dimension, BasisType::GaussLobatto);
	FiniteElementSpace* fes = new FiniteElementSpace(&mesh, fec);
	BilinearForm oldDGTrace(fes);
	BilinearForm maxwellDGTrace(fes);

	std::vector<VectorConstantCoefficient> n = {
	VectorConstantCoefficient(Vector({1.0})),
	};

	const double alpha = 1.0;
	const double beta = 0.0;
	const double gamma = 0.0;

	oldDGTrace.AddInteriorFaceIntegrator(new DGTraceIntegrator(n[0], alpha, beta));
	oldDGTrace.AddBdrFaceIntegrator(new DGTraceIntegrator(n[0], alpha, beta));
	oldDGTrace.Assemble();
	oldDGTrace.Finalize();

	maxwellDGTrace.AddInteriorFaceIntegrator(new MaxwellDGTraceIntegrator(n[0], alpha, beta, gamma));
	maxwellDGTrace.AddBdrFaceIntegrator(new MaxwellDGTraceIntegrator(n[0], alpha, beta, gamma));
	maxwellDGTrace.Assemble();
	maxwellDGTrace.Finalize();

	auto oldDGMatrix = oldDGTrace.SpMat().ToDenseMatrix();
	oldDGMatrix->Print(std::cout);
	std::cout << std::endl;
	auto maxwellDGMatrix = maxwellDGTrace.SpMat().ToDenseMatrix();
	maxwellDGMatrix->Print(std::cout);
	std::cout << std::endl;

	ASSERT_EQ(oldDGMatrix->Height(), maxwellDGMatrix->Width());
	ASSERT_EQ(oldDGMatrix->Height(), maxwellDGMatrix->Height());

	for (int i = 0; i < oldDGMatrix->Width(); i++) {
		for (int j = 0; j < oldDGMatrix->Height(); j++){
			EXPECT_NEAR(maxwellDGMatrix->Elem(i, j), oldDGMatrix->Elem(i, j), tol);
		}
	}
}