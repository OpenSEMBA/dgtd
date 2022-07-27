#include "TestGlobalFunctions.h"


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