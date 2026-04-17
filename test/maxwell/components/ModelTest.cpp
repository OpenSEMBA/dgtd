#include <gtest/gtest.h>

#include "components/Model.h"

using namespace maxwell;
using namespace mfem;

// Helper: 1D mesh, 5 elements, all attribute 1, boundaries attribute 1 (left) and 2 (right).
static Mesh make1DMesh(int n = 5)
{
	return Mesh::MakeCartesian1D(n);
}

// Helper: 1D mesh where half the elements have attribute 2 (for multi-material tests).
static Mesh make1DMeshTwoAttribs(int n = 4)
{
	auto m = Mesh::MakeCartesian1D(n);
	for (int e = n / 2; e < n; e++) {
		m.GetElement(e)->SetAttribute(2);
	}
	m.SetAttributes();
	return m;
}

class ModelTest : public ::testing::Test {
};

// --- numberOfMaterials ---

TEST_F(ModelTest, numberOfMaterials_default)
{
	auto mesh = make1DMesh();
	Model model(mesh);
	EXPECT_EQ(1u, model.numberOfMaterials());
}

TEST_F(ModelTest, numberOfMaterials_singleMaterial)
{
	auto mesh = make1DMesh();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(2.0, 3.0, 0.0));
	Model model(mesh, matInfo);
	EXPECT_EQ(1u, model.numberOfMaterials());
}

TEST_F(ModelTest, numberOfMaterials_twoMaterials)
{
	auto mesh = make1DMeshTwoAttribs();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(2.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(4.0, 9.0, 0.0));
	Model model(mesh, matInfo);
	EXPECT_EQ(2u, model.numberOfMaterials());
}

// --- numberOfBoundaryMaterials ---

TEST_F(ModelTest, numberOfBoundaryMaterials_none)
{
	auto mesh = make1DMesh();
	Model model(mesh);
	EXPECT_EQ(0u, model.numberOfBoundaryMaterials());
}

TEST_F(ModelTest, numberOfBoundaryMaterials_twoBoundaries)
{
	auto mesh = make1DMesh();
	GeomTagToBoundaryInfo bdrInfo;
	bdrInfo.gt2b[1] = BdrCond::PEC;
	bdrInfo.gt2b[2] = BdrCond::PMC;
	Model model(mesh, GeomTagToMaterialInfo{}, bdrInfo);
	EXPECT_EQ(2u, model.numberOfBoundaryMaterials());
}

// --- buildEpsMuPiecewiseVector ---

TEST_F(ModelTest, buildEpsMuPiecewiseVector_E_singleMaterial)
{
	auto mesh = make1DMesh();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(4.0, 9.0, 0.0));
	Model model(mesh, matInfo);

	auto eps = model.buildEpsMuPiecewiseVector(FieldType::E);
	ASSERT_EQ(1, eps.Size());
	EXPECT_DOUBLE_EQ(4.0, eps[0]);
}

TEST_F(ModelTest, buildEpsMuPiecewiseVector_H_singleMaterial)
{
	auto mesh = make1DMesh();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(4.0, 9.0, 0.0));
	Model model(mesh, matInfo);

	auto mu = model.buildEpsMuPiecewiseVector(FieldType::H);
	ASSERT_EQ(1, mu.Size());
	EXPECT_DOUBLE_EQ(9.0, mu[0]);
}

TEST_F(ModelTest, buildEpsMuPiecewiseVector_vacuum)
{
	auto mesh = make1DMesh();
	Model model(mesh);

	auto eps = model.buildEpsMuPiecewiseVector(FieldType::E);
	auto mu  = model.buildEpsMuPiecewiseVector(FieldType::H);
	ASSERT_EQ(1, eps.Size());
	ASSERT_EQ(1, mu.Size());
	EXPECT_DOUBLE_EQ(1.0, eps[0]);
	EXPECT_DOUBLE_EQ(1.0, mu[0]);
}

TEST_F(ModelTest, buildEpsMuPiecewiseVector_twoMaterials)
{
	auto mesh = make1DMeshTwoAttribs();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(2.0, 3.0, 0.0));
	matInfo.gt2m.emplace(2, Material(5.0, 7.0, 0.0));
	Model model(mesh, matInfo);

	auto eps = model.buildEpsMuPiecewiseVector(FieldType::E);
	auto mu  = model.buildEpsMuPiecewiseVector(FieldType::H);

	ASSERT_EQ(2, eps.Size());
	EXPECT_DOUBLE_EQ(2.0, eps[0]); // geomTag 1 -> eps_r 2.0
	EXPECT_DOUBLE_EQ(5.0, eps[1]); // geomTag 2 -> eps_r 5.0

	ASSERT_EQ(2, mu.Size());
	EXPECT_DOUBLE_EQ(3.0, mu[0]); // geomTag 1 -> mu_r 3.0
	EXPECT_DOUBLE_EQ(7.0, mu[1]); // geomTag 2 -> mu_r 7.0
}

// --- buildSigmaPiecewiseVector ---

TEST_F(ModelTest, buildSigmaPiecewiseVector_lossless)
{
	auto mesh = make1DMesh();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(2.0, 1.0, 0.0));
	Model model(mesh, matInfo);

	auto sigma = model.buildSigmaPiecewiseVector();
	ASSERT_EQ(1, sigma.Size());
	EXPECT_DOUBLE_EQ(0.0, sigma[0]);
}

TEST_F(ModelTest, buildSigmaPiecewiseVector_conductive)
{
	auto mesh = make1DMesh();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 50.0));
	Model model(mesh, matInfo);

	auto sigma = model.buildSigmaPiecewiseVector();
	ASSERT_EQ(1, sigma.Size());
	EXPECT_DOUBLE_EQ(50.0, sigma[0]);
}

TEST_F(ModelTest, buildSigmaPiecewiseVector_twoMaterials)
{
	auto mesh = make1DMeshTwoAttribs();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(1.0, 1.0, 25.0));
	Model model(mesh, matInfo);

	auto sigma = model.buildSigmaPiecewiseVector();
	ASSERT_EQ(2, sigma.Size());
	EXPECT_DOUBLE_EQ( 0.0, sigma[0]);
	EXPECT_DOUBLE_EQ(25.0, sigma[1]);
}

// --- Boundary marker assembly ---

TEST_F(ModelTest, boundaryMarkers_PEC_PMC_assembled)
{
	auto mesh = make1DMesh();
	GeomTagToBoundaryInfo bdrInfo;
	bdrInfo.gt2b[1] = BdrCond::PEC;
	bdrInfo.gt2b[2] = BdrCond::PMC;
	Model model(mesh, GeomTagToMaterialInfo{}, bdrInfo);

	auto& bdrMap = model.getBoundaryToMarker();
	EXPECT_EQ(1u, bdrMap.count(BdrCond::PEC));
	EXPECT_EQ(1u, bdrMap.count(BdrCond::PMC));
	EXPECT_EQ(0u, bdrMap.count(BdrCond::SMA));
}

TEST_F(ModelTest, boundaryMarkers_PEC_correctValues)
{
	auto mesh = make1DMesh();
	GeomTagToBoundaryInfo bdrInfo;
	bdrInfo.gt2b[1] = BdrCond::PEC;
	bdrInfo.gt2b[2] = BdrCond::PMC;
	Model model(mesh, GeomTagToMaterialInfo{}, bdrInfo);

	auto& pecMarker = model.getBoundaryToMarker().at(BdrCond::PEC);
	ASSERT_EQ(2, pecMarker.Size());
	EXPECT_EQ(1, pecMarker[0]); // attribute 1 is PEC
	EXPECT_EQ(0, pecMarker[1]); // attribute 2 is not PEC
}

TEST_F(ModelTest, boundaryMarkers_PMC_correctValues)
{
	auto mesh = make1DMesh();
	GeomTagToBoundaryInfo bdrInfo;
	bdrInfo.gt2b[1] = BdrCond::PEC;
	bdrInfo.gt2b[2] = BdrCond::PMC;
	Model model(mesh, GeomTagToMaterialInfo{}, bdrInfo);

	auto& pmcMarker = model.getBoundaryToMarker().at(BdrCond::PMC);
	ASSERT_EQ(2, pmcMarker.Size());
	EXPECT_EQ(0, pmcMarker[0]); // attribute 1 is not PMC
	EXPECT_EQ(1, pmcMarker[1]); // attribute 2 is PMC
}

TEST_F(ModelTest, boundaryMarkers_SMA)
{
	auto mesh = make1DMesh();
	GeomTagToBoundaryInfo bdrInfo;
	bdrInfo.gt2b[1] = BdrCond::SMA;
	Model model(mesh, GeomTagToMaterialInfo{}, bdrInfo);

	auto& bdrMap = model.getBoundaryToMarker();
	EXPECT_EQ(1u, bdrMap.count(BdrCond::SMA));
	EXPECT_EQ(0u, bdrMap.count(BdrCond::PEC));
	EXPECT_EQ(0u, bdrMap.count(BdrCond::PMC));

	auto& smaMarker = bdrMap.at(BdrCond::SMA);
	ASSERT_EQ(2, smaMarker.Size());
	EXPECT_EQ(1, smaMarker[0]);
	EXPECT_EQ(0, smaMarker[1]);
}

TEST_F(ModelTest, invalidGeomTagThrows)
{
	auto mesh = make1DMesh();
	GeomTagToBoundaryInfo bdrInfo;
	bdrInfo.gt2b[0] = BdrCond::PEC; // geomTag 0 is invalid
	EXPECT_ANY_THROW(Model model(mesh, GeomTagToMaterialInfo{}, bdrInfo));
}
