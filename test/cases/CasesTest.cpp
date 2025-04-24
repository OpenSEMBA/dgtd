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
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data, maxwellCase("1D_PEC"), true) };

	GridFunction eOld{ solver.getField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
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

	GridFunction eOld{ solver.getField(E,Y) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 1e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getField(E,Y) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), 1e-2);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
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


TEST_F(CasesTest, 1D_TFSF_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("1D_TFSF"));
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data, maxwellCase("1D_TFSF"), true) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
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
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
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
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
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
		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getFieldProbe(4).getFieldMovies()) {
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
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
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
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
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
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
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
		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getFieldProbe(4).getFieldMovies()) {
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
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
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
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
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
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
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
		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		for (const auto& [t, f] : solver.getFieldProbe(4).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}
}