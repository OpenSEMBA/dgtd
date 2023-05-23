#include <gtest/gtest.h>
#include "maxwell/Solver.h"

#include "SourceFixtures.h"
#include "GlobalFunctions.h"

using Interval = std::pair<double, double>;

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class ProblemDescriptionTest : public ::testing::Test {
protected:
	
};

TEST_F(ProblemDescriptionTest, read_pec_centered_1D)
{
	
}

