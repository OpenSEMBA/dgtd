#include "driver/driver.h"

#include "TestUtils.h"

using json = nlohmann::json;

using namespace mfem;
using namespace maxwell;
using namespace maxwell::driver;

class ExtensiveCasesTest : public ::testing::Test {

};

TEST_F(ExtensiveCasesTest, 2D_PEC_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("2D_PEC"));
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data, maxwellCase("2D_PEC"), true) };

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

	// At the left boundary and right boundaries the electric field should be always close to zero S
	// while the magnetic field should reach +- 2.0 at specific times...
	{
		auto expected_t_half{ 0.5 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_InteriorPEC_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_InteriorBdr_Hesthaven"));
	case_data["model"]["boundaries"][0]["tags"] = { 6, 7 };       // PEC Int
	case_data["model"]["boundaries"][1]["tags"] = { 1, 2, 4, 5 }; // PMC
	case_data["model"]["boundaries"][2]["tags"] = { 3 };          // SMA
	auto solver{ buildSolver(case_data, maxwellCase("2D_InteriorBdr_Hesthaven"), true) };


	GridFunction eOld{ solver.getField(E,Z) };
	auto normOld{ solver.getFields().getNorml2() };

	EXPECT_NE(0.0, normOld);

	double tolerance{ 2e-2 };

	solver.run();

	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	{
		auto expected_t{ 0.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}

	{
		auto expected_t{ 1.5 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	{
		auto expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

	{
		auto expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}
}

TEST_F(ExtensiveCasesTest, 2D_InteriorPMC_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_InteriorBdr_Hesthaven"));
	case_data["model"]["boundaries"][0]["tags"] = { 6, 7 };       // PMC Int
	case_data["model"]["boundaries"][0]["type"] = "PMC";          // PMC Int
	case_data["model"]["boundaries"][1]["tags"] = { 1, 2, 4, 5 }; // PMC
	case_data["model"]["boundaries"][2]["tags"] = { 3 };          // SMA
	case_data["sources"][0]["field_type"] = "H";
	case_data["sources"][0]["polarization"] = { 0.0, 1.0, 0.0 };
	auto solver{ buildSolver(case_data, maxwellCase("2D_InteriorBdr_Hesthaven"), true) };


	GridFunction eOld{ solver.getField(H, Y) };
	auto normOld{ solver.getFields().getNorml2() };

	EXPECT_NE(0.0, normOld);

	double tolerance{ 2e-2 };

	solver.run();

	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	{
		auto expected_t{ 0.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}
	}

	{
		auto expected_t{ 1.5 };
		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}

	{
		auto expected_t{ 4.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

	{
		auto expected_t{ 5.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}
	}
}

TEST_F(ExtensiveCasesTest, 2D_InteriorSMA_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_InteriorBdr_Hesthaven"));
	case_data["model"]["boundaries"][0]["tags"] = { 1, 2, 4, 5 }; // PMC
	case_data["model"]["boundaries"][0]["type"] = "PMC";          // PMC
	case_data["model"]["boundaries"][1]["tags"] = { 3, 6, 7 };    // SMA
	case_data["model"]["boundaries"][1]["type"] = "SMA";          // SMA Int
	case_data["model"]["boundaries"].erase(2);                    // Delete unnecessary entry
	auto solver{ buildSolver(case_data, maxwellCase("2D_InteriorBdr_Hesthaven"), true) };

	GridFunction eOld{ solver.getField(E,Z) };
	auto normOld{ solver.getFields().getNorml2() };

	EXPECT_NE(0.0, normOld);

	double tolerance{ 2e-2 };

	solver.run();

	GridFunction eNew{ solver.getField(E,Z) };
	auto normNew{ solver.getFields().getNorml2() };

	EXPECT_NEAR(0.0, normNew, tolerance);

	{
		auto expected_t{ 0.0 };
		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}

	{
		auto expected_t{ 3.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}

	{
		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}
}

TEST_F(ExtensiveCasesTest, 2D_PEC_Centered_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_PEC"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["order"] = 3;
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_PEC"), true) };

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
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_PEC_Centered_Global)
{
	auto case_data = parseJSONfile(maxwellCase("2D_PEC"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["order"] = 3;
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_PEC"), true) };

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
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}

	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_PEC_Upwind_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_PEC"));
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_PEC"), true) };

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
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}
	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}


TEST_F(ExtensiveCasesTest, 2D_PEC_Upwind_Global)
{
	auto case_data = parseJSONfile(maxwellCase("2D_PEC"));
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_PEC"), true) };

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
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}
	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_PEC_Upwind)
{
	std::string case_name{"2D_PEC"};
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

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
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Ez, tolerance);
			if (std::abs(t - expected_t_half) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}
	}
	//...and at the start and end of the simulation, in the center of the problem,
	// the electric field should be always close to 1.0 due to dimensions of the box.

	for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
		auto expected_t_initial{ 0.0 };
		auto expected_t_final{ 2.0 };
		if (std::abs(t - expected_t_initial) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
		if (std::abs(t - expected_t_final) <= 1e-3) {
			EXPECT_NEAR(1.0, f.Ez, tolerance);
			EXPECT_NEAR(0.0, f.Hy, tolerance);
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Centered_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Centered_Global)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Upwind_TEy_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy"));
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Upwind_TEy_Global)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy"));
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy"));
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Upwind_TEy)
{
	std::string case_name{ "2D_TFSF_TEy" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Upwind_THz)
{
	std::string case_name{ "2D_TFSF_THz" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Centered_TEy_Quads)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy_Quads"));
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy_Quads"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Upwind_TEy_Quads)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy_Quads"));
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy_Quads"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Centered_TEy_Quads_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy_Quads"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy_Quads"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_TFSF_Upwind_TEy_Quads_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_TFSF_TEy_Quads"));
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_TFSF_TEy_Quads"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Hz, tolerance);
				EXPECT_NEAR(2.0, f.Ey, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 2D_Bessel_TEz)
{
	auto case_data = parseJSONfile(maxwellCase("2D_Bessel"));
	auto solver{ buildSolver(case_data, maxwellCase("2D_Bessel"), true) };

	EXPECT_NO_THROW(solver.run());

}

TEST_F(ExtensiveCasesTest, 2D_Bessel_TEz_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("2D_Bessel"));
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_Bessel"), true) };

	EXPECT_NO_THROW(solver.run());

}

TEST_F(ExtensiveCasesTest, 2D_Bessel_TEz_Global)
{
	auto case_data = parseJSONfile(maxwellCase("2D_Bessel"));
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("2D_Bessel"), true) };

	EXPECT_NO_THROW(solver.run());

}

TEST_F(ExtensiveCasesTest, 3D_TFSF_Centered)
{
	auto case_data = parseJSONfile(maxwellCase("3D_TFSF"));
	case_data["solver_options"]["solver_type"] = "centered";
	auto solver{ buildSolver(case_data, maxwellCase("3D_TFSF"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}


TEST_F(ExtensiveCasesTest, 3D_TFSF_Centered_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("3D_TFSF"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("3D_TFSF"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}


TEST_F(ExtensiveCasesTest, 3D_TFSF_Centered_Global)
{
	auto case_data = parseJSONfile(maxwellCase("3D_TFSF"));
	case_data["solver_options"]["solver_type"] = "centered";
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("3D_TFSF"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 3D_TFSF_Upwind)
{
	std::string case_name{ "3D_TFSF" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 3D_TFSF_InteriorPEC_Upwind)
{
	std::string case_name{ "3D_TFSF_InteriorPEC" };
	auto solver{ buildSolverJson(maxwellCase(case_name)) };

	solver.run();

	auto tolerance{ 1e-2 };

	{
		double expected_t{ 1.25 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}


		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}
}


TEST_F(ExtensiveCasesTest, 3D_TFSF_Upwind_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("3D_TFSF"));
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("3D_TFSF"), true) };

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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}
	}

}

TEST_F(ExtensiveCasesTest, 3D_TFSF_InteriorPEC_Upwind_Hesthaven)
{
	auto case_data = parseJSONfile(maxwellCase("3D_TFSF_InteriorPEC"));
	case_data["solver_options"]["hesthaven_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("3D_TFSF_InteriorPEC"), true) };

	solver.run();

	auto tolerance{ 1e-2 };

	{
		double expected_t{ 1.25 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}


		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}
}


TEST_F(ExtensiveCasesTest, 3D_TFSF_InteriorPEC_Upwind_Global)
{
	auto case_data = parseJSONfile(maxwellCase("3D_TFSF_InteriorPEC"));
	case_data["solver_options"]["global_operator"] = true;
	auto solver{ buildSolver(case_data, maxwellCase("3D_TFSF_InteriorPEC"), true) };

	solver.run();

	auto tolerance{ 1e-2 };

	{
		double expected_t{ 1.25 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}


		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(-2.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
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
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(-1.0, f.Ez, tolerance);
				EXPECT_NEAR(-1.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(2).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(3).getFieldMovies()) {
			EXPECT_NEAR(0.0, f.Ex, tolerance);
			EXPECT_NEAR(0.0, f.Ey, tolerance);
			EXPECT_NEAR(0.0, f.Hx, tolerance);
			EXPECT_NEAR(0.0, f.Hz, tolerance);
			if (std::abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}
}