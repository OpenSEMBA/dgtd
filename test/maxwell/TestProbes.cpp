#include "gtest/gtest.h"

#include "maxwell/Probes.h"

using namespace maxwell;
using mfem::DenseMatrix;

class TestProbes : public ::testing::Test {
};

TEST_F(TestProbes, integPointConvChecker_1D)
{
	Probes probes;
	auto pointVec = std::vector<std::vector<double>>({ {0.0}, { 0.5 }, { 1.0 } });
	probes.addPointsProbeToCollection(PointsProbe(E, X, pointVec));
	DenseMatrix testMat(1, 3);
	testMat.Elem(0, 0) = 0.0;
	testMat.Elem(0, 1) = 0.5;
	testMat.Elem(0, 2) = 1.0;
	DenseMatrix probeMat = probes.getPointsProbes().at(0).getIntegPointMat();
	
	for (int i = 0; i < pointVec.at(0).size(); i++) {
		for (int j = 0; j < pointVec.size(); j++) {
			EXPECT_EQ(testMat.Elem(i, j), probeMat.Elem(i, j));
		}
	}
}
TEST_F(TestProbes, integPointDiffDimsVector)
{
	Probes probes;
	auto pointVec = std::vector<std::vector<double>>({ {0.0, 0.5}, {0.5}, {1.0, 0.5, 1.0} });
	ASSERT_ANY_THROW(probes.addPointsProbeToCollection(PointsProbe(E, X, pointVec)));

}

TEST_F(TestProbes, integPointEmptySubvectors)
{
	Probes probes;
	auto pointVec = std::vector<std::vector<double>>({ {},{} });
	ASSERT_ANY_THROW(probes.addPointsProbeToCollection(PointsProbe(E, X, pointVec)));
}

TEST_F(TestProbes, integPointOneSubvector)
{
	Probes probes;
	auto pointVec = std::vector<std::vector<double>>({ {0.5} });
	probes.addPointsProbeToCollection(PointsProbe(E, X, pointVec));
	EXPECT_EQ(1, probes.getPointsProbes().at(0).getIntegPointMat().Size());
	EXPECT_EQ(0.5, probes.getPointsProbes().at(0).getIntegPointMat().Elem(0, 0));
}