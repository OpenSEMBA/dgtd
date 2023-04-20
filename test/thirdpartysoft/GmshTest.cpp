#include "gtest/gtest.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <mfem.hpp>

#include "gtest/gtest.h"

#include "AnalyticalFunctions2D.h"
#include "SourceFixtures.h"
#include "maxwell/Solver.h"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;
using namespace AnalyticalFunctions2D;

class TestGmsh : public ::testing::Test {
protected:
	Probes buildExportProbes()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}
};

TEST_F(TestGmsh, meshDataGmshMSHRead)
{
	/* This mesh includes a physical tag 1 a surface based on a square split in four trianges with a
	point in the center and a physical tag 2 for the four boundary segments.*/
	ASSERT_NO_THROW(Mesh::LoadFromFile("./testData/test.msh", 1, 0));
}

TEST_F(TestGmsh, DISABLED_2DboxwithGmshMesh)
{
	auto mesh = Mesh::LoadFromFile("./testData/test.msh", 1, 0);
	auto fec = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
	auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get(), 1, 0);
	auto model = Model(mesh, AttributeToMaterial{}, AttributeToBoundary{});
	
	maxwell::Solver solver{
		model,
		buildExportProbes(),
		buildGaussianInitialField(E, Z, 0.1, mfem::Vector({0.5,0.5})),
		SolverOptions{}
		.setTimeStep(5e-4)
		.setFinalTime(2.0)
		.setOrder(4)
	};
}

