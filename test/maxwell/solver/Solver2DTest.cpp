#include <gtest/gtest.h>

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"

using namespace maxwell;
using namespace mfem;
using namespace fixtures::sources;

class Solver2DTest : public ::testing::Test
{
protected:
	static const int defaultNumberOfElements_X{3};
	static const int defaultNumberOfElements_Y{3};

	Model buildModel(
		const int nx = defaultNumberOfElements_X,
		const int ny = defaultNumberOfElements_Y,
		const Element::Type elType = Element::Type::TRIANGLE,
		const BdrCond &bdrB = BdrCond::PEC,
		const BdrCond &bdrR = BdrCond::PEC,
		const BdrCond &bdrT = BdrCond::PEC,
		const BdrCond &bdrL = BdrCond::PEC)
	{
		auto msh{Mesh::MakeCartesian2D(nx, ny, elType)};
		return Model(msh, 
			GeomTagToMaterialInfo(), 
			GeomTagToBoundaryInfo(buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL), 
				GeomTagToInteriorBoundary{}));
	}

	Model buildModel(
		const int nx = defaultNumberOfElements_X, const int ny = defaultNumberOfElements_Y,
		const Element::Type elType = Element::Type::TRIANGLE,
		const double sx = 1.0, const double sy = 1.0,
		const BdrCond &bdrB = BdrCond::PEC, const BdrCond &bdrR = BdrCond::PEC,
		const BdrCond &bdrT = BdrCond::PEC, const BdrCond &bdrL = BdrCond::PEC)
	{
		auto msh{Mesh::MakeCartesian2D(nx, ny, elType, false, sx, sy)};
		return Model(msh, 
			GeomTagToMaterialInfo(), 
			GeomTagToBoundaryInfo(buildAttrToBdrMap2D(bdrB, bdrR, bdrT, bdrL), 
				GeomTagToInteriorBoundary{}));
	}

	GeomTagToBoundary buildAttrToBdrMap2D(
		const BdrCond &bdrB, const BdrCond &bdrR,
		const BdrCond &bdrT, const BdrCond &bdrL)
	{
		return {
			{1, bdrB},
			{2, bdrR},
			{3, bdrT},
			{4, bdrL},
		};
	}

	Vector fieldCenter{{0.5, 0.5}};
	
	Probes buildProbes_for_1dot5D()
	{
		// auto probes{buildProbesWithAnExportProbe(20)};
		auto probes{buildProbesEmpty()};
		probes.fieldProbes = {
			FieldProbe{E, Z, {0.0, 0.5}},
			FieldProbe{E, Z, {1.0, 0.5}},
			FieldProbe{H, Y, {0.0, 0.5}},
			FieldProbe{H, Y, {1.0, 0.5}}
		};
		return probes;
	}

	Sources buildPlanewaveForPeriodic()
	{
		return buildPlanewaveInitialField(
			Gaussian{ 0.1, Position({0.0}) },
			Source::Position({0.5, 0.5}),	  // center_
			Source::Polarization(unitVec(Z)), // e polarization_
			Source::Propagation(minusUnitVec(X))	  // propagation direction
		);
	}

	void expectFieldsAreNearAfterEvolution_1dot5D(maxwell::Solver &solver, const double tol=1e-2)
	{
		GridFunction eOld{solver.getField(E, Y)};
		GridFunction hOld{solver.getField(H, Z)};

		solver.run();

		GridFunction eNew{solver.getField(E, Y)};
		GridFunction hNew{solver.getField(H, Z)};

		EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), tol);
		EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), tol);

		// At the left boundary the electric field should be closed to zero and
		// the magnetic field reaches a maximum close to 1.0 or -1.0
		// (the wave splits in two and doubles at the boundary).
		EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tol);
		EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tol);
		EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tol);
		EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tol);
	}

	void expectedFieldsAreNearAfterEvolution_Periodic(maxwell::Solver &solver, const double tol=1e-2)
	{
		GridFunction eOld{solver.getField(E, Y)};
		GridFunction hOld{solver.getField(H, Z)};

		solver.run();

		GridFunction eNew{solver.getField(E, Y)};
		GridFunction hNew{solver.getField(H, Z)};

		EXPECT_NEAR(0.0, eOld.DistanceTo(eNew), tol);
		EXPECT_NEAR(0.0, hOld.DistanceTo(hNew), tol);

		EXPECT_NEAR(1.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tol);
		EXPECT_NEAR(1.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tol);
		EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tol);
		EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tol);
	}
};

TEST_F(Solver2DTest, pec_centered_tris_1dot5D)
{
	// 1dot5D are 2D cases with symmetry along one axis.
	maxwell::Solver solver{
		buildModel(14, 1, Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, mfem::Vector({0.5, 0.5})),
		SolverOptions{}
			.setCentered()
			.setOrder(3)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

TEST_F(Solver2DTest, pec_centered_quads_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			14, 1, Element::Type::QUADRILATERAL, 1.0, 1.0,
			BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setCentered()
			.setOrder(3)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

#ifdef ENABLE_EXTENSIVE_SOLVER_TESTS

TEST_F(Solver2DTest, pec_upwind_tris_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			5, 3, Element::Type::TRIANGLE, 1.0, 1.0,
			BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, Vector({0.5, 0.5}), unitVec(Z)),
		SolverOptions{}
			.setOrder(4)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

TEST_F(Solver2DTest, pec_upwind_quads_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			5, 5, Element::Type::QUADRILATERAL, 1.0, 1.0,
			BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(5e-3) // Automated time estimation fails with quad meshes.
			.setOrder(4)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

TEST_F(Solver2DTest, sma_upwind_tris_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			5, 3, mfem::Element::Type::TRIANGLE, 1.0, 1.0,
			BdrCond::PMC, BdrCond::SMA, BdrCond::PMC, BdrCond::SMA),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setFinalTime(1.0)
			.setOrder(3)
	};

	GridFunction eOld{solver.getField(E, Z)};
	auto zeros{eOld};
	zeros = 0.0;
	EXPECT_TRUE(eOld.DistanceTo(zeros) > 1e-2);

	solver.run();

	double tol{1e-2};
	EXPECT_NEAR(0.0, solver.getField(E,Z).DistanceTo(zeros), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(2).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tol);
}

TEST_F(Solver2DTest, sma_upwind_quads_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			6, 6, mfem::Element::Type::QUADRILATERAL, 1.0, 1.0,
			BdrCond::PMC, BdrCond::SMA, BdrCond::PMC, BdrCond::SMA),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(3e-3)
			.setFinalTime(1.0)
			.setOrder(3)};

	GridFunction eOld{solver.getField(E, Z)};

	auto zeros{eOld};
	zeros = 0.0;
	EXPECT_TRUE(eOld.DistanceTo(zeros) > 1e-2);

	solver.run();

	double tol{1e-3};
	EXPECT_NEAR(0.0, solver.getField(E,Z).DistanceTo(zeros), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(2).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tol);
}

TEST_F(Solver2DTest, periodic_centered_tris)
{
	Mesh m;
	{
		Mesh square{Mesh::MakeCartesian2D(9, 3, Element::TRIANGLE, false, 1.0, 1.0)};
		std::vector<Vector> translations{
			Vector({1.0, 0.0}),
			Vector({0.0, 1.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{m};

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildPlanewaveForPeriodic(),
		SolverOptions{}
			.setTimeStep(5e-3) // Automated time estimation does not work with periodic meshes.
			.setCentered()
			.setOrder(3)
		};

		expectedFieldsAreNearAfterEvolution_Periodic(solver);
}

TEST_F(Solver2DTest, periodic_centered_quads)
{

	Mesh m;
	{
		Mesh square{Mesh::MakeCartesian2D(9, 9, Element::QUADRILATERAL, false, 1.0, 1.0)};
		std::vector<Vector> translations{
			Vector({1.0, 0.0}),
			Vector({0.0, 1.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{m};

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildPlanewaveForPeriodic(),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setCentered()
			.setOrder(3)
	};

	expectedFieldsAreNearAfterEvolution_Periodic(solver);
}

TEST_F(Solver2DTest, periodic_upwind_tris)
{
	Mesh m;
	{
		Mesh square{Mesh::MakeCartesian2D(9, 9, Element::TRIANGLE, false, 1.0, 1.0)};
		std::vector<Vector> translations{
			Vector({1.0, 0.0}),
			Vector({0.0, 1.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{m};

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildPlanewaveForPeriodic(),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setOrder(3)
	};

	expectedFieldsAreNearAfterEvolution_Periodic(solver);
}

TEST_F(Solver2DTest, periodic_upwind_quads)
{
	Mesh m;
	{
		Mesh square{Mesh::MakeCartesian2D(9, 9, Element::QUADRILATERAL, false, 1.0, 1.0)};
		std::vector<Vector> translations{
			Vector({1.0, 0.0}),
			Vector({0.0, 1.0}),
		};
		m = Mesh::MakePeriodic(square, square.CreatePeriodicVertexMapping(translations));
	}

	Model model{m};

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildPlanewaveForPeriodic(),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setOrder(3)
	};

	expectedFieldsAreNearAfterEvolution_Periodic(solver);
}
TEST_F(Solver2DTest, box_resonant_modes_pec_tris)
{
	// Resonant cavity has eigenmodes which have the following frequency.
	// f_{m,n} = \frac{c}{2 \pi \sqrt{\mu_r \varepsilon_r}} 
	//           \sqrt{ (\frac{m \pi}{a})^2 + (\frac{n \pi}{b})^2) }

	maxwell::Solver solver{
		buildModel(
			5, 5, Element::Type::TRIANGLE, 1.0, 1.0,
			BdrCond::PEC, BdrCond::PEC, BdrCond::PEC, BdrCond::PEC),
		buildProbesWithAnExportProbe(10),
		buildResonantModeInitialField(
			E,
			Source::Polarization(unitVec(Z)),
			{1,1}  // m, n modes
		),
		SolverOptions{}
			.setFinalTime(std::sqrt(2)) // Period of resonant cavity is sqrt(2)
			.setOrder(3)
	};


	GridFunction ezOld{solver.getField(E, Z)};
	GridFunction hxOld{solver.getField(H, X)};
	GridFunction hyOld{solver.getField(H, Y)};

	solver.run();

	GridFunction ezNew{solver.getField(E, Z)};
	GridFunction hxNew{solver.getField(H, X)};
	GridFunction hyNew{solver.getField(H, Y)};

	double tol{1e-2};
	EXPECT_NE(0.0, ezOld.DistanceTo(ezNew));
	EXPECT_NE(0.0, hxOld.DistanceTo(hxNew));
	EXPECT_NE(0.0, hyOld.DistanceTo(hyNew));
	EXPECT_NEAR(0.0, ezOld.DistanceTo(ezNew), tol);
	EXPECT_NEAR(0.0, hxOld.DistanceTo(hxNew), tol);
	EXPECT_NEAR(0.0, hyOld.DistanceTo(hyNew), tol);
}

TEST_F(Solver2DTest, interiorPEC_sma_boundaries)
{
	Mesh mesh{Mesh::LoadFromFile(gmshMeshesFolder() + "InteriorPEC2D.msh", 1, 0)};
	mesh.UniformRefinement();
	GeomTagToBoundary attToBdr{{3, BdrCond::PMC}, {4, BdrCond::SMA}};
	GeomTagToInteriorBoundary attToIntBdr{{2, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, attToIntBdr) };

	auto probes{buildProbesWithAnExportProbe(100)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.2, Source::Position({1.0, 0.0}), unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(4.0)
			.setOrder(3)};

	solver.run();
}

TEST_F(Solver2DTest, AutomatedTimeStepEstimator_tri_K2_P3)
{
	auto mesh{Mesh::MakeCartesian2D(1, 1, Element::Type::TRIANGLE)};

	Vector area(mesh.GetNE()), dtscale(mesh.GetNE());
	for (int it = 0; it < mesh.GetNE(); ++it)
	{
		auto el{mesh.GetElement(it)};
		Vector sper(mesh.GetNumFaces());
		sper = 0.0;
		for (int f = 0; f < mesh.GetElement(it)->GetNEdges(); ++f)
		{
			ElementTransformation *T{mesh.GetFaceTransformation(f)};
			const IntegrationRule &ir = IntRules.Get(T->GetGeometryType(), T->OrderJ());
			for (int p = 0; p < ir.GetNPoints(); p++)
			{
				const IntegrationPoint &ip = ir.IntPoint(p);
				sper(it) += ip.weight * T->Weight();
			}
		}
		sper /= 2.0;
		area(it) = mesh.GetElementVolume(it);
		dtscale(it) = area(it) / sper(it);
	}

	auto meshLGL{Mesh::MakeCartesian1D(1, 2.0)};
	DG_FECollection fec{3, 1, BasisType::GaussLegendre};
	FiniteElementSpace fes{&meshLGL, &fec};

	GridFunction nodes(&fes);
	meshLGL.GetNodes(nodes);

	double rmin{abs(nodes(0) - nodes(1))};

	double dt = dtscale.Min() * rmin * 2.0 / 3.0;

	double tol = 1e-8;
	EXPECT_NEAR(0.101761895965867, dt, tol);
}


TEST_F(Solver2DTest, pec_centered_quads_1dot5D_AMR)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.
	The field along Ez is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
	higher than the initial value.*/

	Mesh mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "amr-quad.mesh").c_str(), 1, 0) };
	mesh.UniformRefinement();

	GeomTagToBoundary attToBdr{ {1, BdrCond::PMC}, {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, Vector({0.5, 0.5}), unitVec(Z)),
		SolverOptions{}
			.setCentered()
			.setOrder(3)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver, 1.2e-2);
}

TEST_F(Solver2DTest, pec_quads_1dot5D_AMR)
{
	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.

	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
	A single 2D Gaussian on X and Y is declared along Ez.

	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
	with triangles and squares.
	The field along Ez is extracted before and after the solver calls its run() method and evolves the
	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
	higher than the initial value.*/

	Mesh mesh{ Mesh::LoadFromFile((mfemMeshes2DFolder() + "amr-quad.mesh").c_str(), 1, 0) };
	mesh.UniformRefinement();

	GeomTagToBoundary attToBdr{ {1, BdrCond::PMC}, {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, Vector({0.5, 0.5}), unitVec(Z)),
		SolverOptions{}.setOrder(3)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

#endif