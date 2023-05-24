#pragma once

#include <mfem.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

namespace maxwell {

static Eigen::MatrixXd toEigen(const mfem::DenseMatrix& mat)
{
	Eigen::MatrixXd res(mat.Width(), mat.Height());
	for (int i = 0; i < mat.Width(); i++) {
		for (int j = 0; j < mat.Height(); j++) {
			res(i, j) = mat.Elem(i, j);
		}
	}
	return res;
}

static Eigen::VectorXd toEigenVector(const mfem::Vector& in)
{
	Eigen::VectorXd res;
	res.resize(in.Size());
	for (int i = 0; i < in.Size(); ++i) {
		res(i) = in.Elem(i);
	}
	return res;
}

static Eigen::VectorXcd toEigenComplexVector(const mfem::Vector& in)
{
	Eigen::VectorXcd res;
	res.resize(in.Size());
	for (int i = 0; i < in.Size(); ++i) {
		res(i) = in.Elem(i);
	}
	return res;
}

static mfem::Vector toMFEMVector(const Eigen::VectorXd& in)
{
	mfem::Vector res(int(in.size()));
	for (int i = 0; i < res.Size(); ++i) {
		res(i) = in[i];
	}
	return res;
}

static mfem::SparseMatrix toMFEMSparse(const Eigen::SparseMatrix<double>& sp)
{
	mfem::SparseMatrix res(int(sp.rows()), int(sp.cols()));
	for (int k = 0; k < sp.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(sp, k); it; ++it)
		{
			res.Set(int(it.row()), int(it.col()), it.value());
		}
	res.Finalize();
	return res;
}

}