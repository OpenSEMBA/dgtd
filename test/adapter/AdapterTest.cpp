#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include "ModelAdapter.hpp"

#include <TestUtils.h>

using json = nlohmann::json;

namespace maxwell {

class MaxwellAdapterTest : public ::testing::Test {
};

TEST_F(MaxwellAdapterTest, testFileFound)
{
	EXPECT_NO_THROW(maxwellCase("JSON_Parser_Test"));
}

TEST_F(MaxwellAdapterTest, testFileParsed)
{
	auto file_name{ maxwellCase("JSON_Parser_Test") };
	std::ifstream test_file(file_name);
	EXPECT_NO_THROW(json::parse(test_file));
}

TEST_F(MaxwellAdapterTest, jsonFindsExistingObjects)
{
	auto file_name{ maxwellCase("JSON_Parser_Test") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_TRUE(case_data.contains("solver_options"));
	EXPECT_TRUE(case_data.contains("model"));
	EXPECT_TRUE(case_data.contains("probes"));
	EXPECT_TRUE(case_data.contains("sources"));
}

TEST_F(MaxwellAdapterTest, jsonFindsExistingNestedObjects)
{
	auto file_name{ maxwellCase("JSON_Parser_Test") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_TRUE(case_data["solver_options"].contains("solver_type"));

	EXPECT_TRUE(case_data["model"]["materials"][0].contains("type"));
	EXPECT_TRUE(case_data["model"]["materials"][1].contains("relative_permittivity"));

	EXPECT_TRUE(case_data["probes"].contains("exporter"));
	EXPECT_TRUE(case_data["probes"]["points"][0].contains("position"));

	EXPECT_TRUE(case_data["sources"][0]["magnitude"].contains("delay"));
	EXPECT_TRUE(case_data["sources"][1]["magnitude"].contains("mode"));
}

TEST_F(MaxwellAdapterTest, readsMesh)
{
	auto file_name{ maxwellCase("2D_Parser_BdrAndInterior") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(assembleMeshString(case_data["model"]["filename"]));
	std::string expected{ "./testData/maxwellInputs/2D_Parser_BdrAndInterior/2D_Parser_BdrAndInterior.msh" };
	EXPECT_EQ(expected, assembleMeshString(case_data["model"]["filename"]));
}

TEST_F(MaxwellAdapterTest, adaptsModelObjects)
{
	auto file_name{ maxwellCase("2D_Parser_BdrAndInterior") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(assembleModel(case_data));
	auto model{ assembleModel(case_data) };	

	typedef std::multimap<BdrCond, BoundaryMarker>::iterator MMAPIterator;

	EXPECT_NO_THROW(model.getConstMesh());

	std::pair<MMAPIterator, MMAPIterator> it_pair = 
		model.getBoundaryToMarker().equal_range(BdrCond::PEC);
	auto a{ 0 };
	for (MMAPIterator it = it_pair.first; it != it_pair.second; it++) {
		if (a == 0) {
			mfem::Array<int> exp({ 0,0,0,1,0,0,0 });
			EXPECT_EQ(it->second, exp);
		}
		if (a == 1) {
			mfem::Array<int> exp({ 0,0,0,0,0,1,0 });
			EXPECT_EQ(it->second, exp);
		}
		a++;
	}

	it_pair = model.getBoundaryToMarker().equal_range(BdrCond::PMC);
	a = 0;
	for (MMAPIterator it = it_pair.first; it != it_pair.second; it++) {
		if (a == 0) {
			mfem::Array<int> exp({ 1,0,0,0,0,0,0 });
			EXPECT_EQ(it->second, exp);
		}
		if (a == 1) {
			mfem::Array<int> exp({ 0,0,1,0,0,0,0 });
			EXPECT_EQ(it->second, exp);
		}
		if (a == 2) {
			mfem::Array<int> exp({ 0,0,0,0,1,0,0 });
			EXPECT_EQ(it->second, exp);
		}
		if (a == 3) {
			mfem::Array<int> exp({ 0,0,0,0,0,0,1 });
			EXPECT_EQ(it->second, exp);
		}
		a++;
	}

	it_pair = model.getInteriorBoundaryToMarker().equal_range(BdrCond::PEC);
	for (MMAPIterator it = it_pair.first; it != it_pair.second; it++) {
			mfem::Array<int> exp({ 0,1,0,0,0,0,0 });
			EXPECT_EQ(it->second, exp);
	}
	//EXPECT_EQ(BdrCond::PEC, model.getBoundaryToMarker())
}

}