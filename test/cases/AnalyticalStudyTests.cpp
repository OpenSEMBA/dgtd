#include "CasesFunctions.h"
#include "TestUtils.h"

namespace maxwell{

using namespace mfem;

class AnalyticalStudyTests : public ::testing::Test {

};

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H1_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H1_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P1/2D_Resonant_Box_H1_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H1_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P2/2D_Resonant_Box_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H1_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P3/2D_Resonant_Box_H1_P3.json";

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

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H3_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H1_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P1/2D_Resonant_Box_H3_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H3_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P2/2D_Resonant_Box_H3_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_H3_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_H3_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P3/2D_Resonant_Box_H3_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

// ------------------------------------------------------ //
// ----------- 2D RESONANT BOX LEGENDRE BASIS ----------- //
// ------------------------------------------------------ //

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H1_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H1_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P1/2D_Resonant_Box_BLEG_H1_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H1_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P2/2D_Resonant_Box_BLEG_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H1_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P3/2D_Resonant_Box_BLEG_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H3_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H3_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P1/2D_Resonant_Box_BLEG_H3_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H3_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H3_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P2/2D_Resonant_Box_BLEG_H3_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H3_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H3_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P3/2D_Resonant_Box_BLEG_H3_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

// --------------------------------------------------------- //
// -------------- 3D RESONANT BOX CLOSED BASIS ------------- //
// --------------------------------------------------------- //

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H1_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H1_P1/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H1_P1/3D_Resonant_Box_H1_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H1_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H1_P2/3D_Resonant_Box_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H1_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H1_P3/3D_Resonant_Box_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P4_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

// ------------------------------------------------------ //
// ----------- 3D RESONANT BOX LEGENDRE BASIS ----------- //
// ------------------------------------------------------ //

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H1_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H1_P1/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box_BLEG/3D_Resonant_Box_BLEG_H1_P1/3D_Resonant_Box_BLEG_H1_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H1_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box_BLEG/3D_Resonant_Box_BLEG_H1_P2/3D_Resonant_Box_BLEG_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H1_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box_BLEG/3D_Resonant_Box_BLEG_H1_P3/3D_Resonant_Box_BLEG_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box_BLEG/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box_BLEG/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P4_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/3D_Resonant_Box_BLEG/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}


}