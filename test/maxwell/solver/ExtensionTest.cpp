#include <gtest/gtest.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "solver/SolverExtension.h"
#include "solver/Solver.h"
#include "cases/CasesFunctions.h"
#include "driver/driver.h"
#include "TestUtils.h"


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

TEST_F(SolverExtensionTest, isCorrect_SBC_Properties)
{
    auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SBC_Test"));
    auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SBC_Test"), true);
    for (auto props : global_solver.getSolverOptions().sbc_props){
        if (props.phys_tag == 15){
            ASSERT_EQ(1e-5, props.material_width);
            ASSERT_EQ(15, props.num_of_segments);
            ASSERT_EQ(5, props.order);
            ASSERT_EQ(1e5, props.material.getConductivity());
        }
        else if (props.phys_tag == 18){
            ASSERT_EQ(1e-8, props.material_width);
            ASSERT_EQ(18, props.num_of_segments);
            ASSERT_EQ(8, props.order);
            ASSERT_EQ(1e8, props.material.getConductivity()); // High only so it ends in 8, for 'same-number' easy checks.
        }
    }
}

TEST_F(SolverExtensionTest, targetIds)
{
    size_t order = 2;
    size_t num_of_segments = 10;
    std::vector<NodeId> ids = buildTargetNodeIds(order, num_of_segments);
    ASSERT_EQ(2, ids[0]);
    ASSERT_EQ(3, ids[1]);
    ASSERT_EQ(33, ids[2]);
    ASSERT_EQ(34, ids[3]);
}

TEST_F(SolverExtensionTest, buildTest)
{
    auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SBC_Test"));
    auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SBC_Test"), true);
    for (auto prop : global_solver.getSolverOptions().sbc_props){
        ASSERT_NO_THROW(SBCSolver(global_solver.getModel(), global_solver.getFES(), prop));
    }
}

}