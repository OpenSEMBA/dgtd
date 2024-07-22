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
	auto solver{ buildSolver(case_name) };

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

TEST_F(CasesTest, 2D_PEC_Centered)
{
	auto data{ parseJSONfile("2D_PEC") };
	data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(data) };

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
				EXPECT_NEAR(2.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(2.0, f.Hy, tolerance);
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
	std::string case_data {"2D_PEC"};
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
				EXPECT_NEAR(2.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(2.0, f.Hy, tolerance);
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

TEST_F(CasesTest, 1D_TFSF_Centered)
{
	auto data{ parseJSONfile("1D_TFSF") };
	data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(data) };

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
	std::string case_data{ "1D_TFSF_TEy" };
	auto solver{ buildSolver(case_data) };

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
				EXPECT_NEAR(-2.0, f.Hz, tolerance);
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
				EXPECT_NEAR(-1.0, f.Hz, tolerance);
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
	std::string case_data{ "1D_TFSF_THz" };
	auto solver{ buildSolver(case_data) };

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
				EXPECT_NEAR( 1.0, f.Ey, tolerance);
				EXPECT_NEAR( 1.0, f.Hz, tolerance);
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

TEST_F(CasesTest, 2D_TFSF_Centered)
{
	auto data{ parseJSONfile("2D_TFSF_TEy") };
	data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(data) };

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
	std::string case_data{ "2D_TFSF_TEy" };
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

TEST_F(CasesTest, 2D_TFSF_Upwind_THz)
{
	std::string case_data{ "2D_TFSF_THz" };
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

TEST_F(CasesTest, 2D_NTFF_Box_Upwind)
{
	std::string case_data{ "2D_NTFF_Box" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_SubLambda_Hz)
{
	std::string case_data{ "2D_RCS" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_ThreeLambda_Hz)
{
	std::string case_data{ "2D_RCS_ThreeLambda" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_Salva_Hz)
{
	std::string case_data{ "2D_RCS_Salva" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 2D_RCS_SalvaConfig_Hz)
{
	std::string case_data{ "2D_RCS_SalvaConfig" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 2D_PEC_Bounce45_Upwind)
{
	std::string case_data{ "2D_PEC_Bounce" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 3D_TFSF_Centered)
{
	auto data{ parseJSONfile("3D_TFSF") };
	data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(data) };

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
				EXPECT_NEAR(0.5, f.Hy, tolerance);
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

TEST_F(CasesTest, 3D_TFSF_Upwind)
{
	std::string case_data{ "3D_TFSF" };
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
				EXPECT_NEAR(0.5, f.Hy, tolerance);
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
				EXPECT_NEAR( 1.0, f.Hy, tolerance);
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

TEST_F(CasesTest, 3D_NearToFarField_Upwind)
{
	std::string case_data{ "3D_NearToFarField" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 3D_NearToFarFieldSmaller_Upwind)
{
	std::string case_data{ "3D_NearToFarFieldSmaller" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}

TEST_F(CasesTest, 3D_NTFF_Sphere_Upwind)
{
	std::string case_data{ "3D_RCS" };
	auto solver{ buildSolver(case_data) };

	solver.run();
}
