#include <gtest/gtest.h>

#include <mfem.hpp>

#include <iostream>
#include <filesystem>


namespace maxwell {

using namespace mfem;

class OtherTest : public ::testing::Test {
};

TEST_F(OtherTest, filesystem_path_to_string) 
{
	
	const std::filesystem::path path("testData/mfemMeshes/");

	auto i{ 1 };

	const std::string onedim("testData/mfemMeshes/1D");
	const std::string twodim("testData/mfemMeshes/2D");
	const std::string threedim("testData/mfemMeshes/3D");

	for (auto const& dir_entry : std::filesystem::directory_iterator{ path }) {
	switch (i) {
	case 1:
		EXPECT_EQ(onedim, dir_entry.path().generic_string());
		break;
	case 2:
		EXPECT_EQ(twodim, dir_entry.path().generic_string());
		break;
	case 3:
		EXPECT_EQ(threedim, dir_entry.path().generic_string());
		break;
	}
	i++;
	}

}

TEST_F(OtherTest, filesystem_recursive_iterator_load_mesh) 
{

	const std::filesystem::path path("testData/maxwellInputs/");

	std::filesystem::path target("testData/maxwellInputs/2D_TFSF/2D_TFSF.msh");

	Mesh mesh;

	for (auto const& dir_entry : std::filesystem::recursive_directory_iterator(path)) {
		if (dir_entry == target) {
			mesh = Mesh::LoadFromFile(dir_entry.path().generic_string(), 1, 0);
		}
	}

	EXPECT_EQ(113, mesh.GetNE() + mesh.GetNBE());
	EXPECT_EQ(2, mesh.SpaceDimension());
	EXPECT_EQ(7, mesh.bdr_attributes.Max());

}

TEST_F(OtherTest, filesystem_recursive_iterator_load_data)
{

	const std::filesystem::path base_path("NearToFarFieldExports/middle/");

	std::unique_ptr<GridFunction> Ex, Ey, Ez, Hx, Hy, Hz;
	double time;
	Mesh mesh;

	for (auto const& dir_entry : std::filesystem::recursive_directory_iterator(base_path)) {
		std::ifstream in(dir_entry.path().generic_string());
		if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 4) == "mesh") {
			mesh = Mesh::LoadFromFile(dir_entry.path().generic_string(), 1, 0);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 5) == "Ex.gf") {
			Ex = std::make_unique<GridFunction>(&mesh, in);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 5) == "Ey.gf") {
			Ey = std::make_unique<GridFunction>(&mesh, in);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 5) == "Ez.gf") {
			Ez = std::make_unique<GridFunction>(&mesh, in);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 5) == "Hx.gf") {
			Hx = std::make_unique<GridFunction>(&mesh, in);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 5) == "Hy.gf") {
			Hy = std::make_unique<GridFunction>(&mesh, in);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 5) == "Hz.gf") {
			Hz = std::make_unique<GridFunction>(&mesh, in);
		}
		else if (dir_entry.path().generic_string().substr(dir_entry.path().generic_string().size() - 8) == "time.txt") {
			std::string timeString;
			std::getline(in, timeString);
			time = std::stod(timeString);
		}
		break; //Just for this test.
	}



}


}