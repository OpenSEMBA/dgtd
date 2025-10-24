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
        ASSERT_EQ(1e-3, props.material_width);
        ASSERT_EQ(12, props.num_of_segments);
        ASSERT_EQ(3, props.order);
    }

    const auto& tag2bdr = global_solver.getModel().getGeomTagToIntBoundaryCond();
    for(const auto& [tag, cond] : tag2bdr){
        if (cond == BdrCond::SBC){
            const auto sbc_material = global_solver.getModel().getGeomTagToBoundaryMaterial().at(tag);
            ASSERT_EQ(1e5, sbc_material.getConductivity());
        }
    }


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