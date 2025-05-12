#include <gtest/gtest.h>

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"

#include "TestUtils.h"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class Solver3DTest : public ::testing::Test {
protected:
	static const int defaultNumberOfElements_X{ 3 };
	static const int defaultNumberOfElements_Y{ 3 };
	static const int defaultNumberOfElements_Z{ 3 };
	Model buildModel(
		const int nx = defaultNumberOfElements_X,
		const int ny = defaultNumberOfElements_Y,
		const int nz = defaultNumberOfElements_Z,
		const Element::Type elType = Element::Type::HEXAHEDRON,
		const double sx = 1.0,
		const double sy = 1.0,
		const double sz = 1.0,
		
		const BdrCond& bdr1 = BdrCond::PEC,
		const BdrCond& bdr2 = BdrCond::PEC,
		const BdrCond& bdr3 = BdrCond::PEC,
		const BdrCond& bdr4 = BdrCond::PEC,
		const BdrCond& bdr5 = BdrCond::PEC,
		const BdrCond& bdr6 = BdrCond::PEC) 
	{

		auto msh{ Mesh::MakeCartesian3D(nx, ny, nz, elType, sx, sy, sz) };

		return Model(msh, 
			GeomTagToMaterialInfo{}, 
			GeomTagToBoundaryInfo(buildAttrToBdrMap3D(bdr1, bdr2, bdr3, bdr4, bdr5, bdr6), 
			GeomTagToInteriorBoundary{}));
	}

	static GeomTagToBoundary buildAttrToBdrMap3D(const BdrCond& bdr1, const BdrCond& bdr2, const BdrCond& bdr3, const BdrCond& bdr4, const BdrCond& bdr5, const BdrCond& bdr6)
	{
		return {
			{1, bdr1},
			{2, bdr2},
			{3, bdr3},
			{4, bdr4},
			{5, bdr5},
			{6, bdr6},
		};
	}
};

TEST_F(Solver3DTest, pec_global_1dot5D)
{
	const double tol{ 6e-2 };

	for (const auto& flux : {
				FluxType::Centered, 
				FluxType::Upwind}) {
		for (const auto& elementType : {
					Element::Type::HEXAHEDRON, 
					Element::Type::TETRAHEDRON}) {

			SolverOptions opts;
			opts.setTimeStep(10e-3)
				.setFinalTime(2.0)
				.setOrder(3);		
			opts.evolution.fluxType = flux;
			
			maxwell::Solver solver{
				buildModel(
					10,    1,   1, elementType, 
					1.0, 1.0, 1.0, 
					BdrCond::PEC,BdrCond::PMC,BdrCond::PEC,
					BdrCond::PMC,BdrCond::PEC,BdrCond::PEC
				),
				buildProbesEmpty(),
				buildGaussianInitialField(
					E, 0.1, 
					Source::Position({0.5,0.5,0.5}), 
					unitVec(Z)
				),
				opts
			};

			GridFunction eOld{solver.getField(E, Y)};
			GridFunction hOld{solver.getField(H, Z)};

			solver.run();

			GridFunction eNew{solver.getField(E, Y)};
			GridFunction hNew{solver.getField(H, Z)};

			EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), tol);
			EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), tol);
		}
	}
}

TEST_F(Solver3DTest, periodic_global_cube_hexa)
{ 
	
	const double tol{ 50e-2 };

	for (const auto& flux : {
				FluxType::Centered, 
				FluxType::Upwind}) {
		Mesh m;
		{
			Mesh cube{ Mesh::MakeCartesian3D(6,3,3,Element::HEXAHEDRON,1.0,1.0,1.0) };
			std::vector<Vector> translations{
				Vector({1.0, 0.0, 0.0}),
				Vector({0.0, 1.0, 0.0}),
				Vector({0.0, 0.0, 1.0})
			};
			m = Mesh::MakePeriodic(cube, cube.CreatePeriodicVertexMapping(translations));
		}

		SolverOptions opts;
		opts.setTimeStep(15e-3)
			.setFinalTime(1.0)
			.setOrder(3);
		opts.evolution.fluxType = flux;

		maxwell::Solver solver{
			Model{m},
			buildProbesEmpty(),
			// buildProbesWithAnExportProbe(5),
			buildPlanewaveInitialField(
				Gaussian{0.1, Position({0.0})},
				Source::Position    ({0.5, 0.5, 0.5}), 
				Source::Polarization(unitVec(Y)),
				Source::Propagation(unitVec(X)) 
			),
			opts
		};

		GridFunction eOld{solver.getField(E, Y)};
		GridFunction hOld{solver.getField(H, Z)};

		solver.run();

		GridFunction eNew{solver.getField(E, Y)};
		GridFunction hNew{solver.getField(H, Z)};

		EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), tol);
		EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), tol);
	}
}

TEST_F(Solver3DTest, sma_upwind_hexa_1dot5D)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{{0.0, 0.2, 0.2}},
		PointProbe{{1.0, 0.2, 0.2}},
	};

	maxwell::Solver solver{
		buildModel(
			10, 2, 2, 
			Element::Type::HEXAHEDRON,
			1.0, 0.4, 0.4,
			BdrCond::PEC, BdrCond::PMC, BdrCond::SMA,
			BdrCond::PMC, BdrCond::SMA, BdrCond::PEC),
			probes,
			buildGaussianInitialField(
				E, 0.1,
				Source::Position({ 0.5, 0.2, 0.2 }),
				unitVec(Z)
			),
			SolverOptions{}
			.setTimeStep(7e-3)
			.setFinalTime(1.0)
			.setOrder(3)
	};


	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(0.0, solver.getFields().getNorml2(), tolerance);

	{
		auto expected_t{ 0.5 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.5, f.Ez, tolerance);
				EXPECT_NEAR(0.5, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR( 0.5, f.Ez, tolerance);
				EXPECT_NEAR(-0.5, f.Hy, tolerance);
			}
		}
	}

	{
		auto expected_t{ 1.0 };
		for (const auto& [t, f] : solver.getFieldProbe(0).getFieldMovies()) {
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}

		for (const auto& [t, f] : solver.getFieldProbe(1).getFieldMovies()) {
			if (abs(t - expected_t) <= 1e-3) {
				EXPECT_NEAR(0.0, f.Ez, tolerance);
				EXPECT_NEAR(0.0, f.Hy, tolerance);
			}
		}
	}
}