#include <gtest/gtest.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "solver/SolverExtension.h"
#include "solver/Solver.h"
#include "cases/CasesFunctions.h"
#include "driver/driver.h"
#include "TestUtils.h"
#include "SourceFixtures.h"


namespace maxwell{

using namespace mfem;
using json = nlohmann::json;


class SolverExtensionTest : public ::testing::Test {
protected:

    static json parseJSONfile(const std::string& json_file)
    {
	    std::ifstream test_file(json_file);
	    return json::parse(test_file);
    }

};

TEST_F(SolverExtensionTest, isCorrect_SGBC_Properties)
{
    auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SGBC_Test"));
    auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SGBC_Test"), true);
    double tol = 1e-7;
    for (auto p = 0; p < global_solver.getModel().getSGBCProperties().size(); p++){
        const auto& props = global_solver.getModel().getSGBCProperties()[p];
        ASSERT_EQ(1, props.layers.size());
        const auto& layer = props.layers[0];
        if (p == 0 && props.geom_tags[0] == 15){
            ASSERT_EQ(1e-1, layer.width);
            ASSERT_EQ(5, layer.num_of_segments);
            ASSERT_EQ(1, layer.order);
            ASSERT_NEAR(1e5, layer.material.getConductivity(), tol);
        }
        else if (p == 1 && props.geom_tags[0] == 18){
            ASSERT_EQ(1e-2, layer.width);
            ASSERT_EQ(8, layer.num_of_segments);
            ASSERT_EQ(2, layer.order);
            ASSERT_NEAR(1e8, layer.material.getConductivity(), tol);
        }
    }
}

TEST_F(SolverExtensionTest, buildTest)
{
    auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SGBC_Test"));
    auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SGBC_Test"), true);
    for (const auto prop : global_solver.getModel().getSGBCProperties()){
        ASSERT_NO_THROW(SGBCWrapper::buildSGBCWrapper(prop));
    }
}

}