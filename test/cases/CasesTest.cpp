#include "driver/driver.h"

#include "TestUtils.h"

using json = nlohmann::json;

using namespace mfem;
using namespace maxwell;
using namespace maxwell::driver;

class CasesTest : public ::testing::Test {

};

TEST_F(CasesTest, 1D_PEC_Upwind)
{

	std::string case_name{ "1D_PEC" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

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
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, 
	// the electric field should be always close to one due to symmetries in the problem.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
		if (abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
		}
	}

}

TEST_F(CasesTest, 1D_TFSF_Centered)
{
	auto case_data{ parseJSONfile(maxwellCase("1D_TFSF")) };
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 7.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(2.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(2.0, f.Hy, tolerance);
			}
		}
	}

}


TEST_F(CasesTest, 1D_TFSF_Upwind_TEy)
{
	std::string case_name{ "1D_TFSF_TEy" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ey, tolerance);
				EXPECT_NEAR(2.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ey, tolerance);
				EXPECT_NEAR( 1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 7.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ey, tolerance);
				EXPECT_NEAR(2.0, f.Hz, tolerance);
			}
		}
	}
}

TEST_F(CasesTest, 1D_TFSF_Upwind_THz)
{
	std::string case_name{ "1D_TFSF_THz" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 2e-2 };

	{
		double expected_t{ 7.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ey, tolerance);
				EXPECT_NEAR(2.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ey, tolerance);
				EXPECT_NEAR(1.0, f.Hz, tolerance);
			}
		}
	}

	{
		double expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ey, tolerance);
				EXPECT_NEAR(2.0, f.Hz, tolerance);
			}
		}
	}

}


TEST_F(CasesTest, 2D_PEC_Centered)
{
	auto case_data{ parseJSONfile(maxwellCase("2D_PEC")) };
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data) };

	GridFunction eOld{ solver.getField(E,Z) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 2e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getField(E,Z) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), tolerance);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 1e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(CasesTest, 2D_PEC_Upwind)
{
	std::string case_name{"2D_PEC"};
	auto solver{ buildSolver(maxwellCase(case_name)) };

	GridFunction eOld{ solver.getField(E,Z) };
	auto normOld{ solver.getFields().getNorml2() };

	// Checks fields have been initialized.
	EXPECT_NE(0.0, normOld);

	double tolerance{ 2e-2 };

	solver.run();

	// Checks that field is almost the same as initially because the completion 
	// of a cycle.
	GridFunction eNew{ solver.getField(E,Z) };
	EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), tolerance);

	// Compares all DOFs.
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), 5e-3);

	// At the left boundary and right boundaries the electric field should be always close to zero 
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}
	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(CasesTest, 2D_TFSF_Centered)
{
	auto case_data{ parseJSONfile(maxwellCase("2D_TFSF_TEy")) };
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 1e-2 };

	{
		double expected_t{ 9.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
				EXPECT_NEAR(1.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
				EXPECT_NEAR(1.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(CasesTest, 2D_TFSF_Upwind_TEy)
{
	std::string case_name{ "2D_TFSF_TEy" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 1e-2 };

	{
		double expected_t{ 9.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
				EXPECT_NEAR(1.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
				EXPECT_NEAR(1.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(CasesTest, 2D_TFSF_Upwind_THz)
{
	std::string case_name{ "2D_TFSF_THz" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 1e-2 };

	{
		double expected_t{ 9.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hz, tolerance);
				EXPECT_NEAR(1.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 8.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
				EXPECT_NEAR(1.0, f.Ey, tolerance);
			}
		}
	}

	{
		double expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(CasesTest, 2D_Conductivity_Upwind)
{
	auto case_data{ parseJSONfile(maxwellCase("2D_COND")) };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 2D_Conductivity_Angled_Upwind)
{
	auto case_data{ parseJSONfile(maxwellCase("2D_COND_ANGLED")) };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, DISABLED_2D_NTFF_Box_Upwind) //Unsupported Gmsh element type.
{
	std::string case_name{ "2D_NTFF_Box" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_SubLambda_Hz)
{
	std::string case_name{ "2D_RCS" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_ThreeLambda_Hz)
{
	std::string case_name{ "2D_RCS_ThreeLambda" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_Salva_Hz)
{
	std::string case_name{ "2D_RCS_Salva" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_SalvaConfig_Hz)
{
	std::string case_name{ "2D_RCS_SalvaConfig" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 2D_PEC_Bounce45_Upwind)
{
	std::string case_name{ "2D_PEC_Bounce" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 3D_TFSF_Centered)
{
	auto case_data{ parseJSONfile(maxwellCase("3D_TFSF")) };
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 1e-2 };

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-0.5, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}

TEST_F(CasesTest, 3D_TFSF_Upwind)
{
	std::string case_name{ "3D_TFSF" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	auto normOld{ solver.getFields().getNorml2() };
	EXPECT_EQ(0.0, normOld);

	solver.run();

	double tolerance{ 1e-2 };

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-0.5, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 2.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 6.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}

TEST_F(CasesTest, 3D_NearToFarField_Upwind)
{
	std::string case_name{ "3D_NearToFarField" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 3D_NearToFarFieldSmaller_Upwind)
{
	std::string case_name{ "3D_NearToFarFieldSmaller" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 3D_TFSF_InteriorPEC_Upwind)
{
	std::string case_name{ "3D_TFSF_InteriorPEC" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();

	auto tolerance{ 1e-2 };

	{
		double expected_t{ 1.25 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}


		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

	}

	{
		double expected_t{ 2.25 };

		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}

	{
		double expected_t{ 3.25 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}
}

TEST_F(CasesTest, 3D_NTFF_Sphere_Upwind)
{
	std::string case_name{ "3D_RCS" };
	auto solver{ buildSolver(maxwellCase(case_name)) };

	solver.run();
}

TEST_F(CasesTest, 3D_Conductivity_Upwind)
{
	auto case_data{ parseJSONfile(maxwellCase("3D_COND")) };
	auto solver{ buildSolver(case_data) };

	solver.run();

}

TEST_F(CasesTest, DISABLED_feng_fss)
{
// 	auto probes{ buildProbesWithAnExportProbe(1000) };

// 	std::vector<double> pointR({ 0.01,-0.075,0.06 });
// 	std::vector<double> pointT({ 0.29,-0.075,0.06 });

// 	probes.pointProbes = {
// 		PointProbe{E, X, pointR},
// 		PointProbe{E, X, pointT},		
// 		PointProbe{E, Y, pointR},
// 		PointProbe{E, Y, pointT},		
// 		PointProbe{E, Z, pointR},
// 		PointProbe{E, Z, pointT},
// 		PointProbe{H, X, pointR},
// 		PointProbe{H, X, pointT},
// 		PointProbe{H, Y, pointR},
// 		PointProbe{H, Y, pointT},
// 		PointProbe{H, Z, pointR},
// 		PointProbe{H, Z, pointT}
// 	};

// 	auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "fengfss.msh").c_str(),1,0)};
// 	mesh.Transform(rotateMinus90degAlongZAxis);
// 	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA} };
// 	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

// 	mfem::Vector center_(3);
// 	rotateMinus90degAlongZAxis(Vector({ 0.075,0.075,0.06 }), center_);
// 	mfem::Vector polarization_(3);
// 	rotateMinus90degAlongZAxis(unitVec(Z), polarization_);
	

// 	maxwell::Solver solver{
// 	model,
// 	probes,
// 	buildPlanewaveInitialField(
// 		Gaussian{0.015},
// 		Source::Position({ 0.075,0.075,0.06 }), // center_
// 		Source::Polarization(unitVec(Z)), // e polarization_
// 		Source::Propagation(unitVec(Y)) // propagation direction
// 	),
// 	SolverOptions{}
// 		.setTimeStep(5e-7)
// 		.setFinalTime(0.50)
// 		.setOrder(1)
// 	};

// 	auto normOld{ solver.getFields().getNorml2() };
// 	solver.run();

// 	for (int probeNumber = 0; probeNumber < probes.pointProbes.size(); probeNumber++) {
// 		std::ofstream file("tnf_" + std::to_string(probeNumber) + ".txt"); 
// 		file << "Time and " + std::to_string(probes.pointProbes[probeNumber].getFieldType()) + std::to_string(probes.pointProbes[probeNumber].getDirection()) + "\n";
// 		for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
// 			file << std::to_string(t) + " " + std::to_string(f) + "\n";
// 		}
// 	}

// 	double tolerance{ 1e-2 };
// 	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(CasesTest, DISABLED_feng_fss_symmetry)
{
	// auto probes{ buildProbesWithAnExportProbe(10) };

	// std::vector<double> pointR({ 10.0, 37.5, 30 });
	// std::vector<double> pointT({ 290.0, 37.5, 30 });

	// probes.fieldProbes = {
	// 	FieldProbe{pointR},
	// 	FieldProbe{pointT}
	// };

	// auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "Feng_FSS_Symmetry.msh").c_str(),1,0)};
	// //mesh.Transform(rotateMinus90degAlongZAxis);
	// GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA}};
	// Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	// maxwell::Solver solver{
	// model,
	// probes,
	// buildPlanewaveInitialField(
	// 	Gaussian{16.0},
	// 	Source::Position({ 75.0, 0.0, 0.0 }), // center
	// 	Source::Polarization(unitVec(Z)), // e polarization
	// 	Source::Propagation(unitVec(X)) // propagation direction
	// ),
	// SolverOptions{}
	// 	.setTimeStep(1e-1)
	// 	.setFinalTime(300.0)
	// 	.setOrder(3)
	// };

	// auto normOld{ solver.getFields().getNorml2() };
	// solver.run();

	// for (int probeNumber = 0; probeNumber < probes.fieldProbes.size(); probeNumber++) {
	// 	std::ofstream file(getTestCaseName() + std::to_string(probeNumber) + ".txt");
	// 	file << "Time // Ex // Ey // Ez // Hx // Hy // Hz //""\n";
	// 	for (const auto& fm : solver.getFieldProbe(probeNumber).getFieldMovies()) {
	// 		std::stringstream time, Ex, Ey, Ez, Hx, Hy, Hz;
	// 		time << std::scientific << std::setprecision(7) << (fm.first); 
	// 		Ex << std::scientific << std::setprecision(7) << fm.second.Ex; Ey << std::scientific << std::setprecision(7) << fm.second.Ey; Ez << std::scientific << std::setprecision(7) << fm.second.Ez; 
	// 		Hx << std::scientific << std::setprecision(7) << fm.second.Hx; Hy << std::scientific << std::setprecision(7) << fm.second.Hy; Hz << std::scientific << std::setprecision(7) << fm.second.Hz;
	// 		file << time.str() + " " + Ex.str() + " " + Ey.str() + " " + Ez.str() + " " + Hx.str() + " " + Hy.str() + " " + Hz.str() + "\n";
	// 	}
	// }

	// double tolerance{ 1e-2 };
	// EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(CasesTest, DISABLED_feng_fss_manual)
{
	// auto mesh{ Mesh::LoadFromFile((mfemMeshes3DFolder() + "fengfssmanual.mesh").c_str(),1,0)};
	// Array<Refinement> refinement_list;
	// refinement_list.Append(Refinement(0, 2));
	// refinement_list.Append(Refinement(1, 2));
	// refinement_list.Append(Refinement(2, 2));
	// refinement_list.Append(Refinement(3, 2));
	// refinement_list.Append(Refinement(4, 2));
	// refinement_list.Append(Refinement(5, 2));
	// refinement_list.Append(Refinement(6, 2));
	// refinement_list.Append(Refinement(7, 2));
	// mesh.GeneralRefinement(refinement_list);
	// refinement_list.Append(Refinement(8, 2));
	// refinement_list.Append(Refinement(9, 2));
	// refinement_list.Append(Refinement(10, 2));
	// refinement_list.Append(Refinement(11, 2));
	// refinement_list.Append(Refinement(12, 2));
	// refinement_list.Append(Refinement(13, 2));
	// refinement_list.Append(Refinement(14, 2));
	// refinement_list.Append(Refinement(15, 2));
	// mesh.GeneralRefinement(refinement_list);
	// //refinement_list.Append(Refinement(16, 2));
	// //refinement_list.Append(Refinement(17, 2));
	// //refinement_list.Append(Refinement(18, 2));
	// //refinement_list.Append(Refinement(19, 2));
	// mesh.GeneralRefinement(refinement_list);
	// mesh.Transform(rotateMinus90degAlongZAxis);
	// GeomTagToBoundary attToBdr{ 
	// 	{2, BdrCond::PEC},
	// 	{3, BdrCond::PMC},
	// 	{4, BdrCond::SMA}
	// };
	// GeomTagToInteriorBoundary attToIntBdr{ {5, BdrCond::PEC} };
	// Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, attToIntBdr) };

	// mfem::Vector center(3);
	// rotateMinus90degAlongZAxis(Vector({ 0.15,0.15,0.06 }), center);
	// mfem::Vector polarization(3);
	// rotateMinus90degAlongZAxis(unitVec(Z), polarization);

	// auto probes{ buildProbesWithAnExportProbe(1000) };

	// std::vector<double> pointR({ 0.01,-0.075,0.06 });
	// std::vector<double> pointT({ 0.29,-0.075,0.06 });

	// probes.pointProbes = {
	// 	PointProbe{E, X, pointR},
	// 	PointProbe{E, X, pointT},
	// 	PointProbe{E, Y, pointR},
	// 	PointProbe{E, Y, pointT},
	// 	PointProbe{E, Z, pointR},
	// 	PointProbe{E, Z, pointT},
	// 	PointProbe{H, X, pointR},
	// 	PointProbe{H, X, pointT},
	// 	PointProbe{H, Y, pointR},
	// 	PointProbe{H, Y, pointT},
	// 	PointProbe{H, Z, pointR},
	// 	PointProbe{H, Z, pointT}
	// };

	// maxwell::Solver solver{
	// model,
	// probes,
	// buildPlanewaveInitialField(
	// 	Gaussian{0.015},
	// 	Source::Position({ 0.0 }), // center
	// 	Source::Polarization(unitVec(Z)), // e polarization
	// 	Source::Propagation(unitVec(X)) // propagation direction
	// ),
	// SolverOptions{}
	// 	.setTimeStep(1e-7)
	// 	.setFinalTime(0.0001)
	// 	.setOrder(1)
	// };

	// auto normOld{ solver.getFields().getNorml2() };
	// solver.run();
	
	// double tolerance{ 1e-2 };
	// EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);
}

TEST_F(CasesTest, DISABLED_feng_fss_flat)
{
	// auto probes{ buildProbesWithAnExportProbe(50) };

	// auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_Feng_FSS_Flat.msh").c_str(), 1, 0, true) };
	// GeomTagToBoundary attToBdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::SMA} };
	// GeomTagToInteriorBoundary att2IntCond{ {60, BdrCond::PEC} };
	// Model model(mesh, GeomTagToMaterialInfo{}, GeomTagToBoundaryInfo(attToBdr, att2IntCond));

	// maxwell::Solver solver{
	// model,
	// probes,
	// buildGaussianPlanewave(0.010, 0.1, unitVec(Y), unitVec(X)),
	// SolverOptions{}
	// 	.setTimeStep(9e-4)
	// 	.setFinalTime(1.0)
	// 	.setOrder(3)
	// };

	// solver.run();

}

TEST_F(CasesTest, DISABLED_interiorPEC_fss_hexas)
{
	// auto probes{ buildProbesWithAnExportProbe(2) };

	// std::vector<double> pointR({ 25, 25, 25 });
	// std::vector<double> pointT({ 275, 25, 25 });

	// probes.pointProbes = {
	// 	PointProbe{E, X, pointR},
	// 	PointProbe{E, X, pointT},
	// 	PointProbe{E, Y, pointR},
	// 	PointProbe{E, Y, pointT},
	// 	PointProbe{E, Z, pointR},
	// 	PointProbe{E, Z, pointT},
	// 	PointProbe{H, X, pointR},
	// 	PointProbe{H, X, pointT},
	// 	PointProbe{H, Y, pointR},
	// 	PointProbe{H, Y, pointT},
	// 	PointProbe{H, Z, pointR},
	// 	PointProbe{H, Z, pointT}
	// };

	// auto mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "fsshexas.msh").c_str(),1,0)};
	// GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC},{4,BdrCond::SMA} };
	// Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	// Source::Position center = mfem::Vector({70.0, 0.0, 0.0});

	// maxwell::Solver solver{
	// model,
	// probes,
	// buildPlanewaveInitialField(
	// 	Gaussian{16},
	// 	Source::Position(center), // center
	// 	Source::Polarization(unitVec(Z)), // e polarization
	// 	mfem::Vector(unitVec(X)) // propagation direction
	// ),
	// SolverOptions{}
	// 	.setTimeStep(7.5e-1)
	// 	.setFinalTime(270.0)
	// 	.setOrder(3)
	// };

	// auto normOld{ solver.getFields().getNorml2() };
	// solver.run();

	// for (int probeNumber = 0; probeNumber < probes.pointProbes.size(); probeNumber++) {
	// 	std::ofstream file("fss_sym_" + std::to_string(probeNumber) + ".txt");
	// 	file << "Time and " + std::to_string(probes.pointProbes[probeNumber].getFieldType()) + std::to_string(probes.pointProbes[probeNumber].getDirection()) + "\n";
	// 	for (const auto& [t, f] : solver.getPointProbe(0).getFieldMovie()) {
	// 		file << std::to_string(t) + " " + std::to_string(f) + "\n";
	// 	}
	// }

	// double tolerance{ 1e-2 };
	// EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}
