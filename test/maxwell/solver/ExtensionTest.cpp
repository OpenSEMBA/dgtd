#include <gtest/gtest.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "solver/SolverExtension.h"
#include "solver/Solver.h"
#include "cases/CasesFunctions.h"
#include "driver/driver.h"
#include "TestUtils.h"
#include "SourceFixtures.h"
#include "math/PhysicalConstants.h"


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

    const auto& sgbc_props = global_solver.getModel().getSGBCProperties();
    ASSERT_EQ(2, sgbc_props.size());

    for (const auto& props : sgbc_props) {
        ASSERT_EQ(1, props.geom_tags.size());
        ASSERT_EQ(1, props.layers.size());

        const auto& layer = props.layers[0];

        // Convenience methods must be consistent with single-layer data
        EXPECT_DOUBLE_EQ(layer.width, props.totalWidth());
        EXPECT_EQ(layer.num_of_segments, props.totalSegments());
        EXPECT_EQ(layer.order, props.maxOrder());

        // No sgbc_boundaries defined in the JSON: both sides inactive
        EXPECT_FALSE(props.sgbc_bdr_info.first.isOn);
        EXPECT_FALSE(props.sgbc_bdr_info.second.isOn);

        // Vacuum background (no relative_permittivity/permeability in JSON)
        EXPECT_DOUBLE_EQ(1.0, layer.material.getPermittivity());
        EXPECT_DOUBLE_EQ(1.0, layer.material.getPermeability());

        // n_skin_depths must be positive for conductive materials
        EXPECT_GT(layer.n_skin_depths, 0.0);

        // Tag-specific layer parameters (manual overrides from JSON)
        const size_t tag = props.geom_tags[0];
        if (tag == 15) {
            EXPECT_DOUBLE_EQ(1e-1, layer.width);
            EXPECT_EQ(5, layer.num_of_segments);
            EXPECT_EQ(1, layer.order);
            EXPECT_DOUBLE_EQ(1e5 * physicalConstants::freeSpaceImpedance_SI, layer.material.getConductivity());
        } else if (tag == 18) {
            EXPECT_DOUBLE_EQ(1e-2, layer.width);
            EXPECT_EQ(8, layer.num_of_segments);
            EXPECT_EQ(2, layer.order);
            EXPECT_DOUBLE_EQ(1e8 * physicalConstants::freeSpaceImpedance_SI, layer.material.getConductivity());
        } else {
            FAIL() << "Unexpected SGBC geom_tag: " << tag;
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