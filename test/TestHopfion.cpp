#include "gtest/gtest.h"

#include "Hopfion.h"

#include <fstream>
#include <mfem.hpp>

class TestHopfion : public ::testing::Test {
};

TEST_F(TestHopfion, initialConditionForHopfion) {
	std::string path = "testData/hopfion/";
	
	std::string casename = "hopfion_p1_q1.txt";
	Hopfion hopfion(1, 1);

	std::ifstream inputFile; 
	std::string filename = path + casename;
	inputFile.open(filename);
	if (!inputFile.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}

	while (!inputFile.eof()) {
		double t;
		Hopfion::Vec3 pos, referenceE, referenceH;
		inputFile >> t 
			>> pos[0] >> pos[1] >> pos[2]
			>> referenceE[0] >> referenceE[1] >> referenceE[2] 
			>> referenceH[0] >> referenceH[1] >> referenceH[2];

		Hopfion::FieldEH computed = hopfion.evaluate(t, pos);

		for (std::size_t d = 0; d < 3; d++) {
			EXPECT_NEAR(referenceE[d], computed.first[d], 1e-8);
			EXPECT_NEAR(referenceH[d], computed.second[d], 1e-8);
		}
		
	}
	
}

TEST(Testing, meshCheck) {

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	EXPECT_EQ(nx*ny, mesh.GetNE());

}

TEST(Testing, elementVerticesCheck) {

	int nx = 8; int ny = 8; bool generateEdges = true;
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges);

	std::vector<int> firstElementVerticesVector = { 0, 1, nx+2, nx+1 };
	std::vector<int> lastElementVerticesVector = { nx-1, nx, nx*2+1, nx*2 };
	mfem::Array<int> meshArrayFirstElement;
	mfem::Array<int> meshArrayLastElement;

	mesh.GetElementVertices(0, meshArrayFirstElement);
	mesh.GetElementVertices(nx*ny-1, meshArrayLastElement);

	std::vector<int> vectorFirstElement(meshArrayFirstElement.begin(), meshArrayFirstElement.end());
	std::vector<int> vectorLastElement(meshArrayLastElement.begin(), meshArrayLastElement.end());

	EXPECT_EQ(firstElementVerticesVector, vectorFirstElement);
	EXPECT_EQ(lastElementVerticesVector, vectorLastElement);

}

TEST(Testing, mapElementAndVertex) {

	int nx = 5; int ny = 5; bool generateEdges = true;

	std::vector<int> mapped = mapElementTopLeftVertex(mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, generateEdges));

	EXPECT_EQ(0, mapped[0]);
	EXPECT_EQ(nx * ny - 1, mapped.size() - 1);
	EXPECT_EQ(nx-1,mapped[mapped.size()-1]);
}

