#include "TestUtils.h"

#include "driver/driver.h"

#include <filesystem>
#include <fstream>

using namespace maxwell;
using namespace maxwell::driver;

namespace {

class ScopedCaseDirectory {
public:
	explicit ScopedCaseDirectory(const std::string& case_name)
		: path_(std::filesystem::path(maxwellInputsFolder()) / case_name)
	{
		std::filesystem::remove_all(path_);
		std::filesystem::create_directories(path_);
	}

	~ScopedCaseDirectory()
	{
		std::filesystem::remove_all(path_);
	}

	const std::filesystem::path& path() const
	{
		return path_;
	}

private:
	std::filesystem::path path_;
};

std::string build1DPMLShellMesh(const bool wraps_both_sides)
{
	if (wraps_both_sides) {
		return R"($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
1 1 "PEC_Left"
1 2 "PEC_Right"
$EndPhysicalNames
$Nodes
8
1 0 0 0
2 0.1 0 0
3 0.2 0 0
4 0.4 0 0
5 0.6 0 0
6 0.8 0 0
7 0.9 0 0
8 1.0 0 0
$EndNodes
$Elements
9
1 15 2 1 1 1
2 15 2 2 2 8
3 1 2 2 2 1 2
4 1 2 2 2 2 3
5 1 2 1 1 3 4
6 1 2 1 1 4 5
7 1 2 1 1 5 6
8 1 2 2 2 6 7
9 1 2 2 2 7 8
$EndElements
)";
	}

	return R"($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
1 1 "PEC_Left"
1 2 "PEC_Right"
$EndPhysicalNames
$Nodes
8
1 0 0 0
2 0.1 0 0
3 0.2 0 0
4 0.4 0 0
5 0.6 0 0
6 0.8 0 0
7 0.9 0 0
8 1.0 0 0
$EndNodes
$Elements
9
1 15 2 1 1 1
2 15 2 2 2 8
3 1 2 2 2 1 2
4 1 2 2 2 2 3
5 1 2 1 1 3 4
6 1 2 1 1 4 5
7 1 2 1 1 5 6
8 1 2 1 1 6 7
9 1 2 1 1 7 8
$EndElements
)";
}

std::string build1DCoarsePMLShellMesh()
{
	return R"($MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
1 1 "PEC_Left"
1 2 "PEC_Right"
$EndPhysicalNames
$Nodes
6
1 0 0 0
2 0.2 0 0
3 0.4 0 0
4 0.6 0 0
5 0.8 0 0
6 1.0 0 0
$EndNodes
$Elements
7
1 15 2 1 1 1
2 15 2 2 2 6
3 1 2 2 2 1 2
4 1 2 1 1 2 3
5 1 2 1 1 3 4
6 1 2 1 1 4 5
7 1 2 2 2 5 6
$EndElements
)";
}

json build1DPMLShellCaseData(const std::string& case_name)
{
	return json{
		{ "model", {
			{ "filename", case_name + ".msh" },
			{ "materials", json::array({
				json{{ "tags", json::array({1}) }, { "type", "vacuum" }},
				json{{ "tags", json::array({2}) }, { "type", "PML" }}
			}) },
			{ "boundaries", json::array({
				json{{ "tags", json::array({1, 2}) }, { "type", "PEC" }}
			}) }
		} },
		{ "sources", json::array({
			json{
				{ "type", "initial" },
				{ "field_type", "electric" },
				{ "center", json::array({0.5}) },
				{ "polarization", json::array({0.0, 1.0, 0.0}) },
				{ "dimension", 1 },
				{ "magnitude", json{
					{ "type", "gaussian" },
					{ "spread", 0.1 }
				} }
			}
		}) }
	};
}

void writePMLShellMesh(const ScopedCaseDirectory& case_dir, const std::string& case_name, const bool wraps_both_sides)
{
	std::ofstream mesh_file(case_dir.path() / (case_name + ".msh"));
	mesh_file << build1DPMLShellMesh(wraps_both_sides);
}

void writeCoarsePMLShellMesh(const ScopedCaseDirectory& case_dir, const std::string& case_name)
{
	std::ofstream mesh_file(case_dir.path() / (case_name + ".msh"));
	mesh_file << build1DCoarsePMLShellMesh();
}

}

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

TEST_F(DriverTest, adaptsPMLMaterialObjects)
{
	const std::string case_name = "DriverPMLShellValid";
	ScopedCaseDirectory case_dir(case_name);
	writePMLShellMesh(case_dir, case_name, true);
	auto case_data = build1DPMLShellCaseData(case_name);

	auto model{ buildModel(case_data, case_name, true) };
	auto geom = model.inferPMLBoxGeometry();

	EXPECT_TRUE(model.hasPMLVolumes());
	ASSERT_EQ(1u, model.getPMLRegions().size());
	auto it = model.getPMLRegions().find(2);
	ASSERT_NE(model.getPMLRegions().end(), it);
	EXPECT_TRUE(it->second.matches_vacuum);
	EXPECT_EQ(3u, it->second.grading_order);
	EXPECT_DOUBLE_EQ(1e-6, it->second.target_reflection);
	EXPECT_DOUBLE_EQ(1.0, model.buildEpsMuPiecewiseVector(FieldType::E)[0]);
	EXPECT_DOUBLE_EQ(1.0, model.buildEpsMuPiecewiseVector(FieldType::H)[0]);
	EXPECT_DOUBLE_EQ(0.0, model.buildSigmaPiecewiseVector()[1]);
	EXPECT_NEAR(0.2, geom.inner_min[0], 1e-12);
	EXPECT_NEAR(0.8, geom.inner_max[0], 1e-12);
}

TEST_F(DriverTest, throwsWhenPMLMaterialDoesNotWrapOuterDomain)
{
	const std::string case_name = "DriverPMLShellInvalid";
	ScopedCaseDirectory case_dir(case_name);
	writePMLShellMesh(case_dir, case_name, false);
	auto case_data = build1DPMLShellCaseData(case_name);

	EXPECT_THROW(buildModel(case_data, case_name, true), std::runtime_error);
}

TEST_F(DriverTest, throwsWhenPMLMaterialOverridesVacuumMatch)
{
	auto file_name{ maxwellCase("1D_PEC") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	case_data["model"]["materials"][0]["type"] = "PML";
	case_data["model"]["materials"][0]["relative_permittivity"] = 2.0;

	EXPECT_THROW(buildModel(case_data, file_name, true), std::runtime_error);
}

TEST_F(DriverTest, reportsPMLSetupDiagnostics)
{
	const std::string case_name = "DriverPMLShellDiagnostics";
	ScopedCaseDirectory case_dir(case_name);
	writePMLShellMesh(case_dir, case_name, true);
	auto case_data = build1DPMLShellCaseData(case_name);
	case_data["solver_options"] = json{{ "order", 1 }};

	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(buildModel(case_data, case_name, true));
	const std::string output = testing::internal::GetCapturedStdout();

	EXPECT_NE(output.find("[PML Setup]"), std::string::npos);
	EXPECT_NE(output.find("Normal Cell Size Min"), std::string::npos);
	EXPECT_NE(output.find("Cells Across X (-/+) : 2 / 2"), std::string::npos);
	EXPECT_NE(output.find("Effective DOFs X (-/+) : 4 / 4"), std::string::npos);
	EXPECT_NE(output.find("Sigma Max"), std::string::npos);
}

TEST_F(DriverTest, pmlDiagnosticsDoNotUseFallbackFrequency)
{
	const std::string case_name = "DriverPMLShellNoSpectrum";
	ScopedCaseDirectory case_dir(case_name);
	writePMLShellMesh(case_dir, case_name, true);
	auto case_data = build1DPMLShellCaseData(case_name);
	case_data.erase("sources");

	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(buildModel(case_data, case_name, true));
	const std::string output = testing::internal::GetCapturedStdout();

	EXPECT_NE(output.find("Pulse Max Freq       : unavailable"), std::string::npos);
	EXPECT_EQ(output.find("1 GHz"), std::string::npos);
}

TEST_F(DriverTest, throwsWhenPMLShellIsClearlyUnderresolved)
{
	const std::string case_name = "DriverPMLShellUnderresolved";
	ScopedCaseDirectory case_dir(case_name);
	writeCoarsePMLShellMesh(case_dir, case_name);
	auto case_data = build1DPMLShellCaseData(case_name);
	case_data["solver_options"] = json{{ "order", 1 }};

	EXPECT_THROW(buildModel(case_data, case_name, true), std::runtime_error);
}
