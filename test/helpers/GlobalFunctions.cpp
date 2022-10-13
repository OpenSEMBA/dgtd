#include "GlobalFunctions.h"

#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/Model.h"

using namespace maxwell;
using namespace mfem;
using namespace mfemExtension;

std::unique_ptr<DenseMatrix> toUnique(DenseMatrix* matPtr)
{
	return std::make_unique<DenseMatrix>(*matPtr);
}

Eigen::MatrixXd toEigen(const DenseMatrix& mat)
{
	Eigen::MatrixXd res(mat.Width(), mat.Height());
	for (int i = 0; i < mat.Width(); i++) {
		for (int j = 0; j < mat.Height(); j++) {
			res(i, j) = mat.Elem(i, j);
		}
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

	return toEigen(*toUnique(res.SpMat().ToDenseMatrix()));
}

Eigen::MatrixXd buildInverseMassMatrixEigen(FiniteElementSpace& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(&fes);
	res.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
	res.Assemble();
	res.Finalize();

	return toEigen(*toUnique(res.SpMat().ToDenseMatrix()));
}

Eigen::MatrixXd buildStiffnessMatrixEigen(FiniteElementSpace& fes)
{
	ConstantCoefficient one(1.0);
	BilinearForm res(&fes);
	res.AddDomainIntegrator(new DerivativeIntegrator(one, 0));
	res.Assemble();
	res.Finalize();

	return toEigen(*toUnique(res.SpMat().ToDenseMatrix()));
}

Eigen::MatrixXd	buildNormalPECFluxOperator1D(
	FiniteElementSpace& fes, const std::vector<Direction>& dirVec)
{
	std::vector<Direction> dirs = dirVec;
	AttributeToBoundary attBdr{ {1,BdrCond::PEC},{2,BdrCond::PEC} };

	BilinearForm res(&fes);
	{
		FluxCoefficient c{ 0.0, 1.0 };
		res.AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirs, c.beta));
		//res->AddInteriorFaceIntegrator(new DGTraceIntegrator(*(new VectorConstantCoefficient(Vector(1.0))), 0.0, 0.5));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(fes.GetMesh()->bdr_attributes.Max());
	for (auto const& kv : attBdr) {
		Array<int> bdrMarker(fes.GetMesh()->bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int)kv.first - 1] = 1;

		bdrMarkers[(int)kv.first - 1] = bdrMarker;
		FluxCoefficient c{ 0.0, -2.0 };
		res.AddBdrFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirs, c.beta), bdrMarkers[kv.first - 1]);
		//res->AddBdrFaceIntegrator(new DGTraceIntegrator(*(new VectorConstantCoefficient(Vector(1.0))), 0.0, 2.0), bdrMarkers[kv.first - 1]);
	}
	res.Assemble();
	res.Finalize();

	return toEigen(*toUnique(res.SpMat().ToDenseMatrix()));
}

Eigen::MatrixXd	buildNormalSMAFluxOperator1D(
	FiniteElementSpace& fes, const std::vector<Direction>& dirVec)
{
	std::vector<Direction> dirs = dirVec;
	AttributeToBoundary attBdr{ {1,BdrCond::SMA},{2,BdrCond::SMA} };
	VectorConstantCoefficient one(Vector({ 1.0 }));

	BilinearForm res(&fes);
	{
		//FluxCoefficient c{ 0.0, 1.0 };
		res.AddInteriorFaceIntegrator(new DGTraceIntegrator(one, 0.0, 1.0));
		//res.AddInteriorFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirs, c.beta));
		//res->AddInteriorFaceIntegrator(new DGTraceIntegrator(*(new VectorConstantCoefficient(Vector(1.0))), 0.0, 0.5));
	}

	std::vector<Array<int>> bdrMarkers;
	bdrMarkers.resize(fes.GetMesh()->bdr_attributes.Max());
	for (auto const& kv : attBdr) {
		Array<int> bdrMarker(fes.GetMesh()->bdr_attributes.Max());
		bdrMarker = 0;
		bdrMarker[(int)kv.first - 1] = 1;
		bdrMarkers[(int)kv.first - 1] = bdrMarker;
		//FluxCoefficient c{ 0.0, -1.0 };
		res.AddBdrFaceIntegrator(new DGTraceIntegrator(one, 0.0, -1.0),bdrMarkers[kv.first -1]);
		//res.AddBdrFaceIntegrator(new MaxwellDGTraceJumpIntegrator(dirs, c.beta), bdrMarkers[kv.first - 1]);
		//res->AddBdrFaceIntegrator(new DGTraceIntegrator(*(new VectorConstantCoefficient(Vector(1.0))), 0.0, 2.0), bdrMarkers[kv.first - 1]);
	}
	res.Assemble();
	res.Finalize();

	return toEigen(*toUnique(res.SpMat().ToDenseMatrix()));
}