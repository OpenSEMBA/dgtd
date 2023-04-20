#include "gtest/gtest.h"

#include "SourceFixtures.h"
#include "maxwell/Solver.h"

#include "maxwell/Types.h"
#include "mfem.hpp"
#include "maxwell/Model.h"
#include "maxwell/mfemExtension/BilinearIntegrators.h"
#include "maxwell/mfemExtension/BilinearForm_IBFI.hpp"
#include "maxwell/mfemExtension/LinearIntegrators.h"
#include "maxwell/mfemExtension/LinearForm_IBFI.hpp"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

#ifdef MAXWELL_USE_MPI
class ParSolver2DTest : public ::testing::Test {
protected:
	static const int defaultNumberOfElements_X{ 3 };
	static const int defaultNumberOfElements_Y{ 3 };

	Model buildModel(
		const int nx = defaultNumberOfElements_X,
		const int ny = defaultNumberOfElements_Y,
		const Element::Type elType = Element::Type::TRIANGLE,
		const BdrCond& bdrB = BdrCond::PEC,
		const BdrCond& bdrR = BdrCond::PEC,
		const BdrCond& bdrT = BdrCond::PEC,
		const BdrCond& bdrL = BdrCond::PEC) {
		auto msh{ Mesh::MakeCartesian2D(nx,ny,elType) };
		return Model(msh, AttributeToMaterial{}, buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL));
	}

	AttributeToBoundary buildAttrToBdrMap2D(const BdrCond& bdrB, const BdrCond& bdrR, const BdrCond& bdrT, const BdrCond& bdrL)
	{
		return {
			{1, bdrB},
			{2, bdrR},
			{3, bdrT},
			{4, bdrL},
		};
	}

	Probes buildProbesWithAnExportProbe()
	{
		return { {}, { ExporterProbe{getTestCaseName()} } };
	}

	static std::string getTestCaseName()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}

	Vector fieldCenter{ { 0.5, 0.5 } };
	Source::Polarization zPolarization()
	{
		return Source::Polarization({ 0.0, 0.0, 1.0 });
	}

};

TEST_F(ParSolver2DTest, periodic_tris_mpi)
{
	// Initialize MPI and HYPRE.
	Mpi::Init();
	int num_procs = Mpi::WorldSize();
	int myid = Mpi::WorldRank();
	Hypre::Init();


	auto probes{ buildProbesWithAnExportProbe() };

	probes.exporterProbes[0].visSteps = 1000;

	Mesh m;
	{
		Mesh square{ Mesh::MakeCartesian2D(9, 9, Element::TRIANGLE, false, 2.0, 2.0) };
		std::vector<Vector> translations{
			Vector({2.0, 0.0}),
			Vector({0.0, 2.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	ParMesh parmesh(MPI_COMM_WORLD, m);

	Model model{ parmesh };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{0.2},
			E,
			Source::Position({1.0, 1.0}), // center
			Source::Polarization({0.0, 0.0, 1.0}), // e polarization
			mfem::Vector({1.0, 0.0, 0.0})  // propagation direction
		),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();


}

#endif