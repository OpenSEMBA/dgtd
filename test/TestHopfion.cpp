#include "gtest/gtest.h"

#include "Hopfion.h"

#include <fstream>

class TestHopfion : public ::testing::Test {
};

TEST_F(TestHopfion, initialConditionForHopfion) {
	std::string path = "../../testData/hopfion/";
	
	std::string filename = "hopfion_p1_q1.txt";
	Hopfion hopfion(1, 1);

	std::ifstream inputFile(path + filename);
	while (!inputFile.eof()) {
		double t;
		Hopfion::Vec3 pos, referenceE, referenceH;
		inputFile >> t 
			>> pos[0] >> pos[1] >> pos[2]
			>> referenceE[0] >> referenceE[1] >> referenceE[2] 
			>> referenceH[0] >> referenceH[1] >> referenceH[2];

		Hopfion::FieldEH computed = hopfion.evaluate(t, pos);
		EXPECT_EQ(referenceE, computed.first);
		EXPECT_EQ(referenceH, computed.second);
	}
	
}