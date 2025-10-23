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

    std::unique_ptr<Mesh> buildDefaultSBC1DMesh(size_t num_seg, double sx)
    {
        auto res = Mesh::MakeCartesian1D(num_seg + 2, sx + 2 * (sx/double(num_seg)));
        return std::make_unique<Mesh>(res);
    }

    std::unique_ptr<ParMesh> buildDefaultSBC1DParMesh(Mesh& s_mesh)
    {
        auto res = ParMesh(MPI_COMM_WORLD, s_mesh);
        return std::make_unique<ParMesh>(res);
    }

    std::unique_ptr<DG_FECollection> buildDefaultSBCFEC(size_t order)
    {
        auto res = DG_FECollection(order, 1, BasisType::GaussLobatto);
        return std::make_unique<DG_FECollection>(res);
    }

    std::unique_ptr<ParFiniteElementSpace> buildDefaultSBCParFES(ParMesh& p_mesh, DG_FECollection& fec)
    {
        auto res = ParFiniteElementSpace(&p_mesh, &fec);
        return std::make_unique<ParFiniteElementSpace>(res);
    }

    SBCProperties buildDefaultSBCProperties(size_t num_seg = 10, size_t ord = 1, double width = 1e-4)
    {
        SBCProperties res;
        res.material_width = width;
        res.num_of_segments = num_seg;
        res.order = ord;
        return res;
    }

    static json parseJSONfile(const std::string& json_file)
    {
	    std::ifstream test_file(json_file);
	    return json::parse(test_file);
    }

    // void buildClassDefaultUtilities()
    // {
    //     auto sbcp_ = buildDefaultSBCProperties();
    //     auto s_mesh_ = buildDefaultSBC1DMesh(sbcp_.num_of_segments, sbcp_.material_width);
    //     auto p_mesh_ = buildDefaultSBC1DParMesh(*s_mesh_);
    //     auto fec_ = buildDefaultSBCFEC(sbcp_.order);
    //     auto pfes_ = buildDefaultSBCParFES(*p_mesh_, *fec_);
    // }

    SBCProperties sbcp_;
    std::unique_ptr<Mesh> s_mesh_;
    std::unique_ptr<ParMesh> p_mesh_;
    std::unique_ptr<DG_FECollection> fec_;
    std::unique_ptr<ParFiniteElementSpace> pfes_;

};

TEST_F(SolverExtensionTest, isCorrect_SBC_Properties)
{
    auto case_data = parseJSONfile(maxwellCase("2D_InteriorBoundary_SBC_Test"));
    auto global_solver = driver::buildSolver(case_data, maxwellCase("2D_InteriorBoundary_SBC_Test"), true);
    ASSERT_EQ(1e-3, global_solver.getSolverOptions().sbc_props.material_width);
    ASSERT_EQ(12, global_solver.getSolverOptions().sbc_props.num_of_segments);
    ASSERT_EQ(3, global_solver.getSolverOptions().sbc_props.order);

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
    SBCSolver sbc_solver(global_solver.getModel(), *pfes_, global_solver.getSolverOptions().sbc_props);
}

}