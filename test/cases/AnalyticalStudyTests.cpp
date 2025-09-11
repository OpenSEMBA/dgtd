#include "CasesFunctions.h"
#include "TestUtils.h"

namespace maxwell{

using namespace mfem;

class AnalyticalStudyTests : public ::testing::Test {

};

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P1)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json";

    L2SimDataCalculator calc(data_path, json_file);

}

}