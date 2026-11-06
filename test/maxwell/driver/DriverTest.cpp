#include "TestUtils.h"

#include "driver/driver.h"

#include <filesystem>
#include <fstream>

using namespace maxwell;
using namespace maxwell::driver;

class DriverTest : public ::testing::Test {
};

TEST_F(DriverTest, testFileFound)
{
	EXPECT_NO_THROW(maxwellCase("JSON_Parser_Test"));
}

TEST_F(DriverTest, testFileParsed)
{
	auto file_name{ maxwellCase("JSON_Parser_Test") };
	std::ifstream test_file(file_name);
	EXPECT_NO_THROW(json::parse(test_file));
}

TEST_F(DriverTest, jsonFindsExistingNestedObjects)
{
	auto file_name{ maxwellCase("JSON_Parser_Test") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_TRUE(case_data.contains("solver_options"));
	EXPECT_TRUE(case_data["solver_options"].contains("evolution_operator"));

	EXPECT_TRUE(case_data.contains("model"));
	EXPECT_TRUE(case_data["model"]["materials"][0].contains("type"));
	EXPECT_TRUE(case_data["model"]["materials"][1].contains("relative_permittivity"));

	EXPECT_TRUE(case_data.contains("probes"));
	EXPECT_TRUE(case_data["probes"].contains("exporter"));
	EXPECT_TRUE(case_data["probes"]["point"][0].contains("position"));

	EXPECT_TRUE(case_data.contains("sources"));
	EXPECT_TRUE(case_data["sources"][0]["magnitude"].contains("spread"));
	EXPECT_TRUE(case_data["sources"][1]["magnitude"].contains("mode"));
}

TEST_F(DriverTest, readsMesh)
{
	auto file_name{ maxwellCase("2D_Parser_BdrAndInterior") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(assembleMeshString(case_data["model"]["filename"]));
	std::string expected{ "./testData/maxwellInputs/2D_Parser_BdrAndInterior/2D_Parser_BdrAndInterior.msh" };
	EXPECT_EQ(expected, assembleMeshString(case_data["model"]["filename"]));
}

TEST_F(DriverTest, adaptsModelObjects)
{
	auto file_name{ maxwellCase("2D_Parser_BdrAndInterior") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(buildModel(case_data, file_name, true));
	auto model{ buildModel(case_data, file_name, true) };

	EXPECT_NO_THROW(model.getConstMesh());
	
	// For this specific test problem, we defined the PEC markers on tags 2, 4 and 6. 
	// But tag number 2 will be an interior tag, which will be checked independently.
	// That means our BoundaryToMarker will have tags marked on 4 and 6 for PEC...
	{ 
		auto marker = model.getBoundaryToMarker().find(BdrCond::PEC);
		mfem::Array<int> exp({ 0,0,0,1,0,1,0 });
		EXPECT_EQ(marker->second, exp);
	}
	// ... and 1, 3, 5 and 7 for PMC.
	{
		auto marker = model.getBoundaryToMarker().find(BdrCond::PMC);
		mfem::Array<int> exp({ 1,0,1,0,1,0,1 });
		EXPECT_EQ(marker->second, exp);
	}
	//Whereas our interior boundary PEC is on position 2.
	{
		auto marker = model.getInteriorBoundaryToMarker().find(BdrCond::PEC);
		mfem::Array<int> exp({ 0,1,0,0,0,0,0 });
		EXPECT_EQ(marker->second, exp);
	}
}

TEST_F(DriverTest, adaptsProbeObjects) 
{
	auto file_name{ maxwellCase("1D_PEC") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	// We expect our Adapter will not throw an error while we build the probes...
	EXPECT_NO_THROW(buildProbes(case_data));
	auto probes{ buildProbes(case_data) };

	// ...and as per our problem definition, we expect to find an exporter probe and three field probes.
	EXPECT_EQ(1, probes.exporterProbes.size());
	EXPECT_EQ(3, probes.pointProbes.size());

}

TEST_F(DriverTest, adaptsSourcesObjects)
{
	auto file_name{ maxwellCase("2D_Parser_BdrAndInterior") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(buildSources(case_data));
	auto sources{ buildSources(case_data) };

	EXPECT_EQ(1, sources.size());
}

TEST_F(DriverTest, throwsWhenCaseNamingStyleIsInconsistent)
{
	const std::string case_name = "DriverStyleCheckMismatch";
	const std::filesystem::path case_dir = std::filesystem::path(maxwellInputsFolder()) / case_name;
	const std::filesystem::path json_path = case_dir / (case_name + ".json");

	std::filesystem::create_directories(case_dir);
	std::ofstream test_file(json_path);
	test_file << R"({"model":{"filename":"WrongMeshName.msh"}})";
	test_file.close();

	EXPECT_THROW(buildSolverJson(json_path.string()), std::runtime_error);

	std::filesystem::remove_all(case_dir);
}

namespace {

json loadPecCase()
{
	std::ifstream in(maxwellCase("1D_PEC"));
	return json::parse(in);
}

json debyeMaterial(int tag, double eps_inf, double eps_s, double tau)
{
	return {
		{"tags", json::array({tag})},
		{"debye", {{"eps_inf", eps_inf}, {"eps_s", eps_s}, {"tau", tau}}}
	};
}

} // namespace

TEST_F(DriverTest, debyeRejectsBadMaterials)
{
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"type", "vacuum"}, {"debye", {{"eps_inf", 2.0}, {"eps_s", 4.0}, {"tau", 1e-9}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}}, {"debye", {{"eps_inf", 2.0}, {"eps_s", 4.0}, {"tau", 1e-9}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}}},
			debyeMaterial(1, 2.0, 4.0, 1e-9)
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			debyeMaterial(1, 2.0, 4.0, 1e-9),
			{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		auto mat = debyeMaterial(1, 2.0, 4.0, 1e-9);
		mat["relative_permittivity"] = 3.0;
		case_data["model"]["materials"] = json::array({mat});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({debyeMaterial(1, 0.5, 4.0, 1e-9)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 2.0, 1e-9)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 4.0, 0.0)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"debye", {{"eps_inf", 2.0}, {"eps_s", 4.0}, {"tau", 1e-9}, {"omega_p", 1.0}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["boundaries"][0]["type"] = "SGBC";
		case_data["model"]["boundaries"][0]["material"] = {
			{"debye", {{"eps_inf", 2.0}, {"eps_s", 4.0}, {"tau", 1e-9}}}
		};
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
}

TEST_F(DriverTest, debyeRejectsHesthavenAndSpectral)
{
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 6.0, 1e-9)});
		case_data["solver_options"]["evolution_operator"] = "hesthaven";
		EXPECT_THROW(buildSolver(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 6.0, 1e-9)});
		case_data["solver_options"]["evolution_operator"] = "maxwell";
		EXPECT_THROW(buildSolver(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 6.0, 1e-9)});
		case_data["solver_options"]["spectral"] = true;
		EXPECT_THROW(buildSolver(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
}

TEST_F(DriverTest, debyeStoresEpsInfAndSolverTau)
{
	auto case_data = loadPecCase();
	const double tau_si = 2.0e-9;
	case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 6.0, tau_si)});
	auto model = buildModel(case_data, maxwellCase("1D_PEC"), true);
	ASSERT_TRUE(model.hasDebye());
	const DebyeProperties* pole = model.findDebye(1);
	ASSERT_NE(pole, nullptr);
	EXPECT_DOUBLE_EQ(pole->eps_inf, 2.0);
	EXPECT_DOUBLE_EQ(pole->eps_s, 6.0);
	EXPECT_DOUBLE_EQ(pole->tau_solver, tau_si * physicalConstants::speedOfLight_SI);
	EXPECT_DOUBLE_EQ(model.getGeomTagToMaterial().at(1).getPermittivity(), 2.0);
}

json lorentzMaterial(int tag, double eps_inf, double omega_p, double omega_1, double gamma)
{
	return {
		{"tags", json::array({tag})},
		{"lorentz", {
			{"eps_inf", eps_inf},
			{"omega_p", omega_p},
			{"omega_1", omega_1},
			{"gamma", gamma}
		}}
	};
}

json lorentzMaterialHz(int tag, double eps_inf, double f_p, double f_1, double gamma_hz)
{
	return {
		{"tags", json::array({tag})},
		{"lorentz", {
			{"eps_inf", eps_inf},
			{"f_p", f_p},
			{"f_1", f_1},
			{"gamma", gamma_hz}
		}}
	};
}

TEST_F(DriverTest, lorentzRejectsBadMaterials)
{
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"type", "vacuum"}, {"lorentz", {{"eps_inf", 1.0}, {"omega_p", 1.0}, {"omega_1", 0.0}, {"gamma", 0.0}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}}, {"lorentz", {{"eps_inf", 1.0}, {"omega_p", 1.0}, {"omega_1", 0.0}, {"gamma", 0.0}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		auto both = debyeMaterial(1, 2.0, 4.0, 1e-9);
		both["lorentz"] = {{"eps_inf", 2.0}, {"omega_p", 1.0}, {"omega_1", 0.0}, {"gamma", 0.0}};
		case_data["model"]["materials"] = json::array({both});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			debyeMaterial(1, 2.0, 4.0, 1e-9),
			lorentzMaterial(1, 1.0, 1.0, 0.0, 0.0)
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		auto mat = lorentzMaterial(1, 2.0, 1.0e9, 0.0, 0.0);
		mat["relative_permittivity"] = 3.0;
		case_data["model"]["materials"] = json::array({mat});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterial(1, 0.5, 1.0, 0.0, 0.0)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterial(1, 1.0, 0.0, 0.0, 0.0)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterial(1, 1.0, 1.0, -1.0, 0.0)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterial(1, 1.0, 1.0, 0.0, -1.0)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"lorentz", {{"eps_inf", 1.0}, {"omega_p", 1.0}, {"omega_1", 0.0}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"lorentz", {{"eps_inf", 2.0}, {"f_p", 1.0e8}, {"omega_1", 0.0}, {"gamma", 0.0}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"lorentz", {{"eps_inf", 2.0}, {"f_p", 1.0e8}, {"gamma", 0.0}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterialHz(1, 2.0, 0.0, 1.0e8, 0.0)});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["boundaries"][0]["type"] = "SGBC";
		case_data["model"]["boundaries"][0]["material"] = {
			{"lorentz", {{"eps_inf", 1.0}, {"omega_p", 1.0}, {"omega_1", 0.0}, {"gamma", 0.0}}}
		};
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
}

TEST_F(DriverTest, dispersiveOnPmlIsRejectedAtInit)
{
	auto expectMessage = [](const json& materials) {
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = materials;
		try {
			buildModel(case_data, maxwellCase("1D_PEC"), true);
			FAIL() << "expected dispersive-on-PML rejection";
		} catch (const std::runtime_error& err) {
			EXPECT_STREQ(kDispersiveOnPmlNotAllowed, err.what());
		}
	};
	expectMessage(json::array({
		{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}},
		 {"debye", {{"eps_inf", 2.0}, {"eps_s", 4.0}, {"tau", 1e-9}}}}
	}));
	expectMessage(json::array({
		{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}},
		 {"lorentz", {{"eps_inf", 1.0}, {"omega_p", 1.0}, {"omega_1", 0.0}, {"gamma", 0.0}}}}
	}));
	expectMessage(json::array({
		{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}}},
		debyeMaterial(1, 2.0, 4.0, 1e-9)
	}));
	expectMessage(json::array({
		lorentzMaterial(1, 1.0, 1.0, 0.0, 0.0),
		{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}}}
	}));

	auto case_data = loadPecCase();
	case_data["model"]["materials"] = json::array({debyeMaterial(1, 2.0, 4.0, 1e-9)});
	auto model = buildModel(case_data, maxwellCase("1D_PEC"), true);
	PMLProperties pml;
	pml.geom_tags = {1};
	pml.active_axes = {X};
	try {
		model.setPMLProperties({pml});
		FAIL() << "expected dispersive-on-PML rejection";
	} catch (const std::runtime_error& err) {
		EXPECT_STREQ(kDispersiveOnPmlNotAllowed, err.what());
	}
}

TEST_F(DriverTest, lorentzStoresSolverRates)
{
	auto case_data = loadPecCase();
	const double omega_p = 2.0e9;
	const double omega_1 = 1.0e9;
	const double gamma = 3.0e8;
	case_data["model"]["materials"] = json::array({lorentzMaterial(1, 2.0, omega_p, omega_1, gamma)});
	auto model = buildModel(case_data, maxwellCase("1D_PEC"), true);
	ASSERT_TRUE(model.hasLorentz());
	const LorentzProperties* pole = model.findLorentz(1);
	ASSERT_NE(pole, nullptr);
	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;
	EXPECT_DOUBLE_EQ(pole->eps_inf, 2.0);
	EXPECT_DOUBLE_EQ(pole->omega_p, omega_p * inv_c);
	EXPECT_DOUBLE_EQ(pole->omega_1, omega_1 * inv_c);
	EXPECT_DOUBLE_EQ(pole->gamma, gamma * inv_c);
	EXPECT_DOUBLE_EQ(model.getGeomTagToMaterial().at(1).getPermittivity(), 2.0);
}

TEST_F(DriverTest, gaussianSpreadFromF1eOrSpread)
{
	json spread_only = {{"spread", 0.15}};
	EXPECT_DOUBLE_EQ(assembleGaussianSpread(spread_only, "test"), 0.15);

	const double f_1e = 3.0e8;
	json hz_only = {{"f_1e", f_1e}};
	EXPECT_DOUBLE_EQ(
		assembleGaussianSpread(hz_only, "test"),
		gaussianSpreadForPower1e(f_1e));

	json both = {{"spread", 0.20}, {"f_1e", f_1e}};
	EXPECT_DOUBLE_EQ(
		assembleGaussianSpread(both, "test"),
		gaussianSpreadForPower1e(f_1e));

	json neither = {{"mean", 0.0}};
	EXPECT_THROW(assembleGaussianSpread(neither, "test"), std::runtime_error);
	json bad_f = {{"f_1e", 0.0}};
	EXPECT_THROW(assembleGaussianSpread(bad_f, "test"), std::runtime_error);
}

TEST_F(DriverTest, lorentzHzStyleMatchesRadStyle)
{
	auto case_data = loadPecCase();
	const double f_p = 2.5e8;
	const double f_1 = 2.0e8;
	const double gamma_hz = 2.5e7;
	case_data["model"]["materials"] = json::array({lorentzMaterialHz(1, 2.0, f_p, f_1, gamma_hz)});
	auto model = buildModel(case_data, maxwellCase("1D_PEC"), true);
	ASSERT_TRUE(model.hasLorentz());
	const LorentzProperties* pole = model.findLorentz(1);
	ASSERT_NE(pole, nullptr);
	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;
	const double two_pi = 2.0 * M_PI;
	EXPECT_DOUBLE_EQ(pole->eps_inf, 2.0);
	EXPECT_DOUBLE_EQ(pole->omega_p, two_pi * f_p * inv_c);
	EXPECT_DOUBLE_EQ(pole->omega_1, two_pi * f_1 * inv_c);
	EXPECT_DOUBLE_EQ(pole->gamma, two_pi * gamma_hz * inv_c);
}

TEST_F(DriverTest, lorentzRejectsHesthavenAndSpectral)
{
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterial(1, 1.0, 1.0e9, 0.0, 0.0)});
		case_data["solver_options"]["evolution_operator"] = "hesthaven";
		EXPECT_THROW(buildSolver(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({lorentzMaterial(1, 1.0, 1.0e9, 0.0, 0.0)});
		case_data["solver_options"]["spectral"] = true;
		EXPECT_THROW(buildSolver(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
}

TEST_F(DriverTest, drudeAndSpatialProfilesParseCorrectly)
{
	const double inv_c = 1.0 / physicalConstants::speedOfLight_SI;

	// 1. Drude with RAM-C 61 km stagnation-line profiles stretched over [0, 1] on 1D_PEC
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{
				{"tags", {1}},
				{"drude", {
					{"eps_inf", 1.0},
					{"omega_p", {{"type", "ramc_61km_omega_p"}}},
					{"nu_c", {{"type", "ramc_61km_nu_c_gupta"}}}
				}}
			}
		});
		auto model = buildModel(case_data, maxwellCase("1D_PEC"), true);
		ASSERT_TRUE(model.hasLorentz());
		const LorentzProperties* pole = model.findLorentz(1);
		ASSERT_NE(pole, nullptr);
		EXPECT_TRUE(pole->hasSpatialProfile());
		EXPECT_TRUE(pole->hasRegionBounds());
		EXPECT_DOUBLE_EQ(pole->region_min_x, 0.0);
		EXPECT_DOUBLE_EQ(pole->region_max_x, 1.0);
		EXPECT_DOUBLE_EQ(pole->omega_1, 0.0);

		mfem::Vector pt_mesh(1);
		pt_mesh(0) = 0.5;
		mfem::Vector pt_ramc(1);
		pt_ramc(0) = -0.012 * (1.0 - pt_mesh(0));
		EXPECT_NEAR(pole->evalOmegaP(pt_mesh), ramc61kmStagnationOmegaP(pt_ramc) * inv_c, 1e-12);
		EXPECT_NEAR(pole->evalGamma(pt_mesh), ramc61kmStagnationNuC_Gupta(pt_ramc) * inv_c, 1e-12);
	}

	// 2. Liu closure preset + Gaussian and Gaussian-derivative profiles auto-mapped via cutoff_factor
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{
				{"tags", {1}},
				{"lorentz", {
					{"eps_inf", 1.5},
					{"omega_p", {
						{"type", "gaussian"},
						{"peak", 4.0e9},
						{"cutoff_factor", 1.0e-2},
						{"offset", 1.0e8}
					}},
					{"omega_1", {
						{"type", "positive_gaussian_derivative"},
						{"peak", 2.0e9},
						{"cutoff_factor", 1.0e-2},
						{"offset", 5.0e8}
					}},
					{"nu_c", {{"type", "ramc_61km_nu_c_liu"}}}
				}}
			}
		});
		auto model = buildModel(case_data, maxwellCase("1D_PEC"), true);
		const LorentzProperties* pole = model.findLorentz(1);
		ASSERT_NE(pole, nullptr);
		EXPECT_TRUE(pole->hasSpatialProfile());

		// At endpoints x = 0.0 and x = 1.0 (xi = 0 and xi = 1), both profiles decay to offset + cutoff_factor * peak
		mfem::Vector left_pt(1), right_pt(1), center_pt(1);
		left_pt(0) = 0.0;
		center_pt(0) = 0.5;
		right_pt(0) = 1.0;
		EXPECT_NEAR(pole->evalOmegaP(left_pt), (1.0e8 + 1.0e-2 * 4.0e9) * inv_c, 1e-12);
		EXPECT_NEAR(pole->evalOmegaP(center_pt), (1.0e8 + 4.0e9) * inv_c, 1e-12);
		EXPECT_NEAR(pole->evalOmegaP(right_pt), (1.0e8 + 1.0e-2 * 4.0e9) * inv_c, 1e-12);

		EXPECT_NEAR(pole->evalOmega1(left_pt), (5.0e8 + 1.0e-2 * 2.0e9) * inv_c, 1e-12);
		EXPECT_NEAR(pole->evalOmega1(center_pt), 5.0e8 * inv_c, 1e-12);
		EXPECT_NEAR(pole->evalOmega1(right_pt), (5.0e8 + 1.0e-2 * 2.0e9) * inv_c, 1e-12);

		// Check that both lobes of positive_gaussian_derivative (at z = -1 and z = +1) reach offset + peak
		const auto* deriv_prof = dynamic_cast<const PositiveGaussianDerivativeProfile*>(pole->omega_1_profile.get());
		ASSERT_NE(deriv_prof, nullptr);
		const double z_cut = deriv_prof->zCut();
		mfem::Vector lobe1_pt(1), lobe2_pt(1);
		lobe1_pt(0) = 0.5 * (1.0 - 1.0 / z_cut);
		lobe2_pt(0) = 0.5 * (1.0 + 1.0 / z_cut);
		EXPECT_NEAR(pole->evalOmega1(lobe1_pt), (5.0e8 + 2.0e9) * inv_c, 1e-12);
		EXPECT_NEAR(pole->evalOmega1(lobe2_pt), (5.0e8 + 2.0e9) * inv_c, 1e-12);

		mfem::Vector pt_mesh(1);
		pt_mesh(0) = 0.75;
		mfem::Vector pt_ramc(1);
		pt_ramc(0) = -0.012 * (1.0 - pt_mesh(0));
		EXPECT_NEAR(pole->evalGamma(pt_mesh), ramc61kmStagnationNuC_Liu(pt_ramc) * inv_c, 1e-12);
	}

	// 3. Reject invalid combinations on drude
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}}, {"type", "PML"}, {"active_axes", {"X"}},
			 {"drude", {{"omega_p", 1.0e9}, {"nu_c", 1.0e7}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}
	{
		auto case_data = loadPecCase();
		case_data["model"]["materials"] = json::array({
			{{"tags", {1}},
			 {"lorentz", {{"eps_inf", 1.0}, {"omega_p", 1.0e9}, {"omega_1", 0.0}, {"gamma", 0.0}}},
			 {"drude", {{"omega_p", 1.0e9}, {"nu_c", 1.0e7}}}}
		});
		EXPECT_THROW(buildModel(case_data, maxwellCase("1D_PEC"), true), std::runtime_error);
	}

	// 4. Parse full 1D_Drude_RAMC_61km case JSON
	{
		const std::string case_path = maxwellCase("1D_Drude_RAMC_61km");
		auto case_data = parseJSONfile(case_path);
		auto model = buildModel(case_data, case_path, false);
		EXPECT_TRUE(model.hasLorentz());
		const LorentzProperties* pole = model.findLorentz(3);
		ASSERT_NE(pole, nullptr);
		EXPECT_TRUE(pole->hasSpatialProfile());
	}
}
