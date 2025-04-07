#include "TestUtils.h"

#include "driver/driver.h"

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
	EXPECT_TRUE(case_data["solver_options"].contains("solver_type"));

	EXPECT_TRUE(case_data.contains("model"));
	EXPECT_TRUE(case_data["model"]["materials"][0].contains("type"));
	EXPECT_TRUE(case_data["model"]["materials"][1].contains("relative_permittivity"));

	EXPECT_TRUE(case_data.contains("probes"));
	EXPECT_TRUE(case_data["probes"].contains("exporter"));
	EXPECT_TRUE(case_data["probes"]["point"][0].contains("position"));

	EXPECT_TRUE(case_data.contains("sources"));
	EXPECT_TRUE(case_data["sources"][0]["magnitude"].contains("delay"));
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
