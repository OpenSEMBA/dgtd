#include <gtest/gtest.h>

#include "math/EigenMfemTools.h"

using namespace maxwell;

class EigenMfemToolsTest : public ::testing::Test {
};

TEST_F(EigenMfemToolsTest, vectorRoundTrip)
{
	mfem::Vector mfemV({1.0, 2.5, -3.7, 0.0, 4.2});

	auto eigenV = toEigenVector(mfemV);
	auto back = toMFEMVector(eigenV);

	ASSERT_EQ(mfemV.Size(), back.Size());
	for (int i = 0; i < mfemV.Size(); i++) {
		EXPECT_DOUBLE_EQ(mfemV[i], back[i]);
	}
}

TEST_F(EigenMfemToolsTest, vectorConversion_sizes)
{
	mfem::Vector empty;
	auto eigenEmpty = toEigenVector(empty);
	EXPECT_EQ(0, eigenEmpty.size());

	mfem::Vector single({42.0});
	auto eigenSingle = toEigenVector(single);
	ASSERT_EQ(1, eigenSingle.size());
	EXPECT_DOUBLE_EQ(42.0, eigenSingle(0));
}

TEST_F(EigenMfemToolsTest, denseMatrixConversion)
{
	mfem::DenseMatrix mat(2, 3);
	mat(0, 0) = 1.0; mat(0, 1) = 2.0; mat(0, 2) = 3.0;
	mat(1, 0) = 4.0; mat(1, 1) = 5.0; mat(1, 2) = 6.0;

	auto eigenMat = toEigen(mat);

	ASSERT_EQ(2, eigenMat.rows());
	ASSERT_EQ(3, eigenMat.cols());
	EXPECT_DOUBLE_EQ(1.0, eigenMat(0, 0));
	EXPECT_DOUBLE_EQ(2.0, eigenMat(0, 1));
	EXPECT_DOUBLE_EQ(3.0, eigenMat(0, 2));
	EXPECT_DOUBLE_EQ(4.0, eigenMat(1, 0));
	EXPECT_DOUBLE_EQ(5.0, eigenMat(1, 1));
	EXPECT_DOUBLE_EQ(6.0, eigenMat(1, 2));
}

TEST_F(EigenMfemToolsTest, complexVectorConversion)
{
	mfem::Vector mfemV({1.0, -2.5, 3.0});

	auto eigenCV = toEigenComplexVector(mfemV);

	ASSERT_EQ(3, eigenCV.size());
	for (int i = 0; i < mfemV.Size(); i++) {
		EXPECT_DOUBLE_EQ(mfemV[i], eigenCV(i).real());
		EXPECT_DOUBLE_EQ(0.0, eigenCV(i).imag());
	}
}

TEST_F(EigenMfemToolsTest, sparseMatrixConversion)
{
	Eigen::SparseMatrix<double> eigenSp(3, 3);
	eigenSp.insert(0, 0) = 1.0;
	eigenSp.insert(1, 2) = 5.0;
	eigenSp.insert(2, 1) = -3.0;
	eigenSp.makeCompressed();

	auto mfemSp = toMFEMSparse(eigenSp);

	ASSERT_EQ(3, mfemSp.NumRows());
	ASSERT_EQ(3, mfemSp.NumCols());
	EXPECT_EQ(3, mfemSp.NumNonZeroElems());

	// Verify nonzero entries via GetRow
	mfem::Array<int> cols;
	mfem::Vector vals;

	mfemSp.GetRow(0, cols, vals);
	ASSERT_EQ(1, cols.Size());
	EXPECT_EQ(0, cols[0]);
	EXPECT_DOUBLE_EQ(1.0, vals[0]);

	mfemSp.GetRow(1, cols, vals);
	ASSERT_EQ(1, cols.Size());
	EXPECT_EQ(2, cols[0]);
	EXPECT_DOUBLE_EQ(5.0, vals[0]);

	mfemSp.GetRow(2, cols, vals);
	ASSERT_EQ(1, cols.Size());
	EXPECT_EQ(1, cols[0]);
	EXPECT_DOUBLE_EQ(-3.0, vals[0]);
}

TEST_F(EigenMfemToolsTest, sparseMatrixConversion_identity)
{
	int n = 4;
	Eigen::SparseMatrix<double> identity(n, n);
	identity.setIdentity();

	auto mfemId = toMFEMSparse(identity);

	ASSERT_EQ(n, mfemId.NumRows());
	ASSERT_EQ(n, mfemId.NumCols());
	EXPECT_EQ(n, mfemId.NumNonZeroElems());

	// Verify diagonal entries
	mfem::Array<int> cols;
	mfem::Vector vals;
	for (int i = 0; i < n; i++) {
		mfemId.GetRow(i, cols, vals);
		ASSERT_EQ(1, cols.Size());
		EXPECT_EQ(i, cols[0]);
		EXPECT_DOUBLE_EQ(1.0, vals[0]);
	}
}
