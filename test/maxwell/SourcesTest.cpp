#include "gtest/gtest.h"

#include "Solver.h"
#include "Utils.h"
#include "SourceFixtures.h"
#include "GlobalFunctions.h"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

using Solver = maxwell::Solver;

class SourcesTest : public ::testing::Test {
protected:
	
};

TEST_F(SourcesTest, planewave)
{
	Planewave pw{Gaussian{ 0.2, mfem::Vector({-1.5}) }, unitVec(Y), unitVec(X)};

	const double TOL{ 1e-6 };

	// Follows maximum.
	EXPECT_NEAR(1.0, pw.eval(Source::Position({ -1.5 }), 0.0, E, Y), TOL);
	EXPECT_NEAR(1.0, pw.eval(Source::Position({  0.0 }), 1.5, E, Y), TOL);
	
	// H has same magnitude as E.
	for (auto x{0.0}; x < 3.0; x += 0.1) {
		EXPECT_NEAR(
			pw.eval(Source::Position({ x }), 0.0, E, Y), 
			pw.eval(Source::Position({ x }), 0.0, H, Z), 
			TOL
		);
	}	

}
