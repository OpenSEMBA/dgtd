#include "gtest/gtest.h"

#include "maxwell/Probes.h"

using namespace maxwell;

class TestMaxwellProbes : public ::testing::Test {
};

TEST_F(TestMaxwellProbes, integPointConvChecker_1D)
{
	Probes probes;
	auto pointVec = std::vector<std::vector<double>>({ {0.0}, { 0.5 }, { 1.0 } });
	probes.addProbeToVector(Probe(E, X, pointVec));
	DenseMatrix testMat(1, 3);
	testMat.Elem(0, 0) = 0.0;
	testMat.Elem(0, 1) = 0.5;
	testMat.Elem(0, 2) = 1.0;
	DenseMatrix probeMat = probes.getProbeVector().at(0).getIntegPointMat();
	
	for (int i = 0; i < pointVec.at(0).size(); i++) {
		for (int j = 0; j < pointVec.size(); j++) {
			EXPECT_EQ(testMat.Elem(i, j), probeMat.Elem(i, j));
		}
	}
}

TEST_F(TestMaxwellProbes, integPointFailVector)
{
	Probes probes;
	auto pointVec = std::vector<std::vector<double>>({ {0.0, 0.5}, {0.5}, {1.0, 0.5, 1.0} });
	probes.addProbeToVector(Probe(E, X, pointVec));

	//Expect throw test, should fail, diff dimms on vector.
}