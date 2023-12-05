#include <gtest/gtest.h>

#include <fftw3.h>

#include <mfem.hpp>
#include <components/Types.h>

#include <iostream>
#include <filesystem>

#include <math.h>
#include <complex>


namespace maxwell {

using namespace mfem;

std::string getGridFunctionString(const FieldType& f, const Direction& d) 
{
	switch (f) {
	case E:
		switch (d) {
		case X:
			return "Ex.gf";
		case Y:
			return "Ey.gf";
		case Z:
			return "Ez.gf";
		}
	case H:
		switch (d) {
		case X:
			return "Hx.gf";
		case Y:
			return "Hy.gf";
		case Z:
			return "Hz.gf";
		}
	}
}

std::string getGridFunctionPathForType(const std::string& path, const FieldType& f, const Direction& d)
{
	return path + "/" + getGridFunctionString(f, d);
}

Mesh getRCSMesh(const std::string& path)
{
	std::ifstream in(path + "/mesh");
	return Mesh(in);
}

class RCSTest : public ::testing::Test {
public:
	GridFunction getGridFunction(Mesh& m, const std::string& path, const FieldType& f, const Direction& d) 
	{
		auto filepath{ getGridFunctionPathForType(path, f, d) };
		std::ifstream in(filepath);
		if (!in) {
			throw std::runtime_error("File could not be opened in readGridFunctionFromFile, verify path.");
		}
		GridFunction gf(&m, in);
		return gf;
	}

};

TEST_F(RCSTest, readMeshFromFile)
{
	EXPECT_NO_THROW(getRCSMesh("testData/rcsInputs"));
}

TEST_F(RCSTest, readGridFunctionsFromFile)
{
	std::string path{ "testData/rcsInputs" };
	auto m{ getRCSMesh(path) };
	EXPECT_NO_THROW(getGridFunction(m, path, E, X));
}

TEST_F(RCSTest, explorerTest)
{
	std::string path{ "testData/rcsInputs" };
	auto m{ getRCSMesh(path) };
	auto gf{ getGridFunction(m, path, E, X) };
}
}