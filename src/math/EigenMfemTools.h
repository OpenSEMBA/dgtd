#pragma once

#include <mfem.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues> 

namespace maxwell {

inline Eigen::MatrixXd toEigen(const mfem::DenseMatrix& mat)
{
	Eigen::MatrixXd res(mat.Height(), mat.Width());
	for (int i = 0; i < mat.Height(); i++) {
		for (int j = 0; j < mat.Width(); j++) {
			res(i, j) = mat.Elem(i, j);
		}
	}
	return res;
}

inline Eigen::VectorXd toEigenVector(const mfem::Vector& in)
{
	Eigen::VectorXd res;
	res.resize(in.Size());
	for (int i = 0; i < in.Size(); ++i) {
		res(i) = in.Elem(i);
	}
	return res;
}

inline Eigen::VectorXcd toEigenComplexVector(const mfem::Vector& in)
{
	Eigen::VectorXcd res;
	res.resize(in.Size());
	for (int i = 0; i < in.Size(); ++i) {
		res(i) = in.Elem(i);
	}
	return res;
}

inline mfem::Vector toMFEMVector(const Eigen::VectorXd& in)
{
	mfem::Vector res(int(in.size()));
	for (int i = 0; i < res.Size(); ++i) {
		res(i) = in[i];
	}
	return res;
}

inline mfem::SparseMatrix toMFEMSparse(const Eigen::SparseMatrix<double>& sp)
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