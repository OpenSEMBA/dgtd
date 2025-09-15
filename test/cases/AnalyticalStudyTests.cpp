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

// ------------------------------------------------------ //
// ----------- 2D RESONANT BOX LEGENDRE BASIS ----------- //
// ------------------------------------------------------ //

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P1_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_BLEG_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P2_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_BLEG_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P3_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Box_BLEG_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P1_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Box_BLEG_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P1_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Box_BLEG_H2_P1/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P2_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Box_BLEG_H2_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Box_BLEG_H2_P3_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Box_BLEG_H2_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

// ------------------------------------------------------ //
// ----------- 2D RESONANT CIRCLE CLOSED BASIS ----------- //
// ------------------------------------------------------ //

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P2_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Circle_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P2/2D_Resonant_Circle_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P3_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Circle_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P3/2D_Resonant_Circle_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P4_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Circle_H1_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P4/2D_Resonant_Circle_H1_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Circle_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P2/2D_Resonant_Circle_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Circle_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P3/2D_Resonant_Circle_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P4_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Circle_H1_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P4/2D_Resonant_Circle_H1_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P2_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Circle_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P2/2D_Resonant_Circle_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P3_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Circle_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P3/2D_Resonant_Circle_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_H1_P4_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Circle_H1_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_H1_P4/2D_Resonant_Circle_H1_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

// --------------------------------------------------------- //
// ----------- 2D RESONANT CIRCLE LEGENDRE BASIS ----------- //
// --------------------------------------------------------- //

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P2_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Circle_BLEG_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P2/2D_Resonant_Circle_BLEG_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P3_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Circle_BLEG_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P3/2D_Resonant_Circle_BLEG_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P4_cuda)
{
    std::string data_path = "./Exports/cuda-8/2D_Resonant_Circle_BLEG_H1_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P4/2D_Resonant_Circle_BLEG_H1_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Circle_BLEG_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P2/2D_Resonant_Circle_BLEG_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Circle_BLEG_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P3/2D_Resonant_Circle_BLEG_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P4_mpi)
{
    std::string data_path = "./Exports/mpi-8/2D_Resonant_Circle_BLEG_H1_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P4/2D_Resonant_Circle_BLEG_H1_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P2_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Circle_BLEG_H1_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P2/2D_Resonant_Circle_BLEG_H1_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P3_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Circle_BLEG_H1_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P3/2D_Resonant_Circle_BLEG_H1_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 2D_Resonant_Circle_BLEG_H1_P4_single)
{
    std::string data_path = "./Exports/single-core/2D_Resonant_Circle_BLEG_H1_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Circle_BLEG_H1_P4/2D_Resonant_Circle_BLEG_H1_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

// --------------------------------------------------------- //
// -------------- 3D RESONANT BOX CLOSED BASIS ------------- //
// --------------------------------------------------------- //

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P2_cuda)
{
    std::string data_path = "./Exports/cuda-8/3D_Resonant_Box_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P3_cuda)
{
    std::string data_path = "./Exports/cuda-8/3D_Resonant_Box_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P4_cuda)
{
    std::string data_path = "./Exports/cuda-8/3D_Resonant_Box_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P4_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P2_single)
{
    std::string data_path = "./Exports/single-core/3D_Resonant_Box_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P3_single)
{
    std::string data_path = "./Exports/single-core/3D_Resonant_Box_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_H0_P4_single)
{
    std::string data_path = "./Exports/single-core/3D_Resonant_Box_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

// ------------------------------------------------------ //
// ----------- 3D RESONANT BOX LEGENDRE BASIS ----------- //
// ------------------------------------------------------ //

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P2_cuda)
{
    std::string data_path = "./Exports/cuda-8/3D_Resonant_Box_BLEG_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P3_cuda)
{
    std::string data_path = "./Exports/cuda-8/3D_Resonant_Box_BLEG_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

// This is here for completitude but this case does not fit in the current GPU due to the open nature of the basis and the size that comes from it.
// TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P4_cuda)
// {
//     std::string data_path = "./Exports/cuda-8/3D_Resonant_Box_BLEG_H0_P4/";
//     std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json";

//     RMSDataCalculator calc(data_path, json_file);

// }

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P2_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P3_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P4_mpi)
{
    std::string data_path = "./Exports/mpi-8/3D_Resonant_Box_BLEG_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P2_single)
{
    std::string data_path = "./Exports/single-core/3D_Resonant_Box_BLEG_H0_P2/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P3_single)
{
    std::string data_path = "./Exports/single-core/3D_Resonant_Box_BLEG_H0_P3/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json";

    RMSDataCalculator calc(data_path, json_file);

}

TEST_F(AnalyticalStudyTests, 3D_Resonant_Box_BLEG_H0_P4_single)
{
    std::string data_path = "./Exports/single-core/3D_Resonant_Box_BLEG_H0_P4/";
    std::string json_file = "./testData/maxwellInputs/2D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json";

    RMSDataCalculator calc(data_path, json_file);

}


}