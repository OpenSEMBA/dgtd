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
