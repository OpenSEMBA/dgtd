#include "CasesFunctions.h"
#include "TestUtils.h"

namespace maxwell{

using namespace mfem;

class AnalyticalStudyTests : public ::testing::Test {

};

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P1_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P2_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P2/2D_Resonant_Box_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P3_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P3/2D_Resonant_Box_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P2/2D_Resonant_Box_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P3/2D_Resonant_Box_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P1_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Box_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P2_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Box_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P2/2D_Resonant_Box_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H2_P3_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Box_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P3/2D_Resonant_Box_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

}