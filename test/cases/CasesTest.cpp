#include "driver/driver.h"

#include "TestUtils.h"

using json = nlohmann::json;

using namespace mfem;
using namespace maxwell;
using namespace maxwell::driver;

class CasesTest : public ::testing::Test {

};

TEST_F(CasesTest, 1D_PEC_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("1D_PEC"));
    case_data["solver_options"]["upwind_alpha"] = 0.0;
	auto solver{ buildSolver(case_data, maxwellCase("1D_PEC"), true) };

	GridFunction eOld{ solver.getConstField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getConstField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}

}

TEST_F(CasesTest, 1D_PEC_Upwind)
{

	std::string case_name{ "1D_PEC" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

	GridFunction eOld{ solver.getConstField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getConstField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}

}


TEST_F(CasesTest, 1D_SMA_Upwind)
{
	auto case_data = parseJSONfile(maxwellCase("1D_PEC"));
    case_data["model"]["boundaries"][0]["type"] = "SMA";
	case_data["solver_options"]["final_time"] = 1.0;
	auto solver{ buildSolver(case_data, maxwellCase("1D_PEC"), true) };

	GridFunction eOld{ solver.getConstField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getConstField(E,Y) };
	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), 1e-2);

}


TEST_F(CasesTest, 1D_MultiCondition_Upwind)
{
	SGBCProperties sbcp;
	SGBCLayer layer(Material(20.0, 1.0, 20.0), 2.0);
	layer.order = 4;
	layer.num_of_segments = 10;
	sbcp.layers.push_back(layer);
	SGBCBoundaries sbcp_bdrs;
	sbcp_bdrs.second.isOn = true;
	sbcp_bdrs.second.bdrCond = BdrCond::PEC;
	sbcp.sgbc_bdr_info = sbcp_bdrs;

	auto total_segments = sbcp.totalSegments();
	auto total_width = sbcp.totalWidth();
	auto mesh = mfem::Mesh::MakeCartesian1D(total_segments + 2, total_width + 2 * total_width / total_segments);
    mesh.AddBdrPoint(1, 3);
    mesh.AddBdrPoint(mesh.GetNV() - 2, 4);
    mesh.SetAttribute(0, 2);
    mesh.SetAttribute(mesh.GetNE() - 1, 2);
    mesh.bdr_attributes = mfem::Array<int>({1, 2, 3, 4}); 
    mesh.FinalizeMesh();
    int* partitioning = mesh.GeneratePartitioning(1);

    Material vacuum = buildVacuumMaterial();
    GeomTagToMaterial geom_tag_sgbc_mat{{1, vacuum}, {2, vacuum}};
    GeomTagToInteriorBoundary gt2ib;
    if (sbcp_bdrs.first.isOn){
        gt2ib[3] = sbcp_bdrs.first.bdrCond;
    }
    if (sbcp_bdrs.second.isOn){
        gt2ib[4] = sbcp_bdrs.second.bdrCond;
    }
	GeomTagToBoundary gt2b;
    gt2b[1] = BdrCond::SMA;
    gt2b[2] = BdrCond::SMA;
    GeomTagToBoundaryInfo gtbdr(gt2b, gt2ib);
    Model model = Model(mesh, GeomTagToMaterialInfo(geom_tag_sgbc_mat, GeomTagToBoundaryMaterial{}), gtbdr, partitioning);

	Probes probes;
    probes.exporterProbes.resize(1);
    ExporterProbe ep;
    ep.name = "1D_MultiCondition_Upwind";
    ep.visSteps = 100;
    probes.exporterProbes.at(0) = ep;
    Sources sources;
	mfem::Vector gaussianCenter(1);
	gaussianCenter = 0.5;
	double spread = 0.1;
	Gaussian gauss{ spread, gaussianCenter, 1 };
	sources.add(std::make_unique<InitialField>(gauss, E, mfem::Vector({0.0,1.0,0.0}), gaussianCenter));
	sources.add(std::make_unique<InitialField>(gauss, H, mfem::Vector({0.0,0.0,1.0}), gaussianCenter));
    SolverOptions opts;
    opts.setOrder(sbcp.maxOrder());
    opts.setUpwindAlpha(1.0);
    opts.setODEType(ode_type::RK4); 
	opts.setFinalTime(10.0);
	opts.setTimeStep(1e-3);
    opts.setExportEO(true);
    std::cout << "Assembling SGBC Solvers: " << std::endl;

    auto solver = maxwell::Solver(model, probes, sources, opts);

	solver.run();

}


TEST_F(CasesTest, 1D_TFSF_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("1D_TFSF"));
    case_data["solver_options"]["upwind_alpha"] = 0.0;
	auto solver{ buildSolver(case_data, maxwellCase("1D_TFSF"), true) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 4.2 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 10.2 };
		for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 7.2 };
		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getPointProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getPointProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}

TEST_F(CasesTest, 1D_TFSF_Upwind)
{
	std::string case_name{ "1D_TFSF" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 4.2 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 10.2 };
		for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 7.2 };
		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getPointProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getPointProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}

TEST_F(CasesTest, 1D_TFSF_Upwind_Negative_K)
{
	auto case_data = parseJSONfile(maxwellCase("1D_TFSF"));
	case_data["probes"]["point"][0]["position"] = {2.99999};
	case_data["probes"]["point"][2]["position"] = { -2.99999 };
	case_data["sources"][0]["propagation"] = { -1.0, 0.0, 0.0 };
	case_data["sources"][0]["magnitude"]["mean"] = { 5.0 };
	auto solver{ buildSolver(case_data, maxwellCase("1D_TFSF"), true) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getPointProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getPointProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getPointProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getPointProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}