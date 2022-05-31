#include "gtest/gtest.h"

#include "maxwell/Probes.h"

using namespace maxwell;

class TestMaxwellProbes : public ::testing::Test {
};

TEST_F(TestMaxwellProbes, integPointConvChecker_1D)
{
	Probes probe;
	auto pointVec = std::vector<std::vector<double>>({ {0.0}, { 0.5 }, { 1.0 } });
	probe.addProbeToVector(Probe(E, X, pointVec));
	DenseMatrix testMat(1, 3);
	testMat.Elem(0, 0) = 0.0;
	testMat.Elem(0, 1) = 0.5;
	testMat.Elem(0, 2) = 1.0;
	DenseMatrix probeMat = probe.getProbeVector().at(0).getIntegPointMat();
	
	for (int i = 0; i < pointVec.at(0).size(); i++) {
		for (int j = 0; j < pointVec.size(); j++) {
			EXPECT_EQ(testMat.Elem(i, j), probeMat.Elem(i, j));
		}
	}

}