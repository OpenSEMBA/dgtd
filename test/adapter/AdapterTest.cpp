#pragma once
#include <fstream>
#include <nlohmann/json.hpp>

#include "MaxwellAdapter.hpp"
#include "ModelAdapter.hpp"
#include "ProbesAdapter.hpp"
#include "SourcesAdapter.hpp"
#include "SolverOptsAdapter.hpp"

#include <TestUtils.h>

using json = nlohmann::json;
using namespace mfem;

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
	EXPECT_TRUE(case_data["probes"]["field"][0].contains("position"));

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

	EXPECT_NO_THROW(buildModel(case_data));
	auto model{ buildModel(case_data) };	

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

TEST_F(MaxwellAdapterTest, adaptsProbeObjects) 
{
	auto file_name{ maxwellCase("1D_PEC_Centered") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(buildProbes(case_data));
	auto probes{ buildProbes(case_data) };

	EXPECT_EQ(1, probes.exporterProbes.size());
	EXPECT_EQ(1, probes.fieldProbes.size());

}

TEST_F(MaxwellAdapterTest, adaptsSourcesObjects)
{
	auto file_name{ maxwellCase("2D_Parser_BdrAndInterior") };
	std::ifstream test_file(file_name);
	auto case_data = json::parse(test_file);

	EXPECT_NO_THROW(buildSources(case_data));
	auto sources{ buildSources(case_data) };

	EXPECT_EQ(1, sources.size());
}

TEST_F(MaxwellAdapterTest, 1D_PEC_Centered) 
{
	std::string case_name{ "1D_PEC_Centered" };
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

	// At the left boundary the electric field should be always close to zero...
	for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
		EXPECT_NEAR(0.0, f.Ey, tolerance);
	}

}

TEST_F(MaxwellAdapterTest, 1D_PEC_Upwind)
{

	std::string case_name{ "1D_PEC_Upwind" };
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

	// At the left boundary the electric field should be always close to zero...
	for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
		EXPECT_NEAR(0.0, f.Ey, tolerance);
	}

}

TEST_F(MaxwellAdapterTest, 2D_PEC_Centered)
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

	// At the left boundary the electric field should be always close to zero...
	for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
		EXPECT_NEAR(0.0, f.Ex, tolerance);
		EXPECT_NEAR(0.0, f.Ey, tolerance);
		EXPECT_NEAR(0.0, f.Ez, tolerance);
	}

}

TEST_F(MaxwellAdapterTest, 2D_PEC_Upwind)
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

	// At the left boundary the electric field should be always close to zero...
	for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
		EXPECT_NEAR(0.0, f.Ex, tolerance);
		EXPECT_NEAR(0.0, f.Ey, tolerance);
		EXPECT_NEAR(0.0, f.Ez, tolerance);
	}

}


}
