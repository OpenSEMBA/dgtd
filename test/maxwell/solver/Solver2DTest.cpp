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
		auto probes{buildProbesWithAnExportProbe(20)};
		// auto probes{buildProbesEmpty()};
		probes.pointProbes = {
			PointProbe{E, Z, {0.0, 0.5}},
			PointProbe{E, Z, {1.0, 0.5}},
			PointProbe{H, Y, {0.0, 0.5}},
			PointProbe{H, Y, {1.0, 0.5}}
		};
		return probes;
	}

	Sources buildPlanewaveForPeriodic()
	{
		return buildPlanewaveInitialField(
			Gaussian{0.1},
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

class Solver2DSpectralTest : public Solver2DTest
{
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

TEST_F(Solver2DSpectralTest, DISABLED_pec_centered_tris_1dot5D_spectral)
{

	auto probes{buildProbesWithAnExportProbe()};
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}};

	maxwell::Solver solver{
		buildModel(14, 1, Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(3)
			.setSpectralEO(true)};

	auto normOld{solver.getFields().getNorml2()};
	solver.run();

	double tolerance{1e-2};
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
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

TEST_F(Solver2DSpectralTest, DISABLED_pec_centered_quads_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}};

	maxwell::Solver solver{
		buildModel(14, 1, Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(3)
			.setSpectralEO()};

	auto normOld{solver.getFields().getNorml2()};
	solver.run();

	double tolerance{1e-2};
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
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

	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "amr-quad.mesh").c_str(), 1, 0)};
	mesh.UniformRefinement();

	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, Vector({0.5, 0.5}), unitVec(Z)),
		SolverOptions{}
			.setCentered()
			.setOrder(3)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver,1.2e-2);
}

TEST_F(Solver2DSpectralTest, DISABLED_pec_tris_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}};

	maxwell::Solver solver{
		buildModel(14, 1, Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setOrder(3)
			.setSpectralEO()
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
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

	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "amr-quad.mesh").c_str(), 1, 0)};
	mesh.UniformRefinement();

	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	maxwell::Solver solver{
		model,
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, Vector({0.5, 0.5}), unitVec(Z)),
		SolverOptions{}.setOrder(3)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

TEST_F(Solver2DTest, pec_tris_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			5, 3, Element::Type::TRIANGLE, 1.0, 1.0,
			BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, Vector({0.5, 0.5}), unitVec(Z)),
		SolverOptions{}.setOrder(5)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

TEST_F(Solver2DTest, pec_quads_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			5, 5, Element::Type::QUADRILATERAL, 1.0, 1.0,
			BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		buildProbes_for_1dot5D(),
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(5e-3) // Automated time estimation fails with quad meshes.
			.setOrder(5)
	};

	expectFieldsAreNearAfterEvolution_1dot5D(solver);
}

TEST_F(Solver2DSpectralTest, DISABLED_pec_quads_1dot5D_spectral)
{

	Probes probes;
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0, 0.5}}};

	maxwell::Solver solver{
		buildModel(14, 1, Element::Type::QUADRILATERAL, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(2.0)
			.setOrder(3)
			.setSpectralEO()};

	auto normOld{solver.getFields().getNorml2()};
	solver.run();

	double tolerance{1e-2};
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	// At the left boundary the electric field should be closed to zero and
	// the magnetic field reaches a maximum close to 1.0 or -1.0
	// (the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}

TEST_F(Solver2DTest, sma_tris_1dot5D)
{
	maxwell::Solver solver{
		buildModel(
			10, 3, mfem::Element::Type::TRIANGLE, 1.0, 1.0,
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

	double tol{1e-3};
	EXPECT_NEAR(0.0, solver.getField(E,Z).DistanceTo(zeros), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(2).findFrameWithMin().second), tol);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tol);
}

TEST_F(Solver2DTest, sma_quads_1dot5D)
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

TEST_F(Solver2DSpectralTest, DISABLED_periodic_centered_tris_spectral_and_base_comparison)
{

	Probes probes;

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
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(3)};

	maxwell::Solver solverSpectral{
		model,
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(3)
			.setSpectralEO()};

	Vector zeroVec{solver.getField(E, Z).Size()};
	zeroVec = 0.0;
	double tolerance{1e-5};
	for (int i = 0; i < zeroVec.Size(); ++i)
	{
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E, Z).Elem(i) - solverSpectral.getField(E, Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, Y).Elem(i) - solverSpectral.getField(H, Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, X).Elem(i) - solverSpectral.getField(H, X).Elem(i), tolerance);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < zeroVec.Size(); ++i)
	{
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E, Z).Elem(i) - solverSpectral.getField(E, Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, Y).Elem(i) - solverSpectral.getField(H, Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, X).Elem(i) - solverSpectral.getField(H, X).Elem(i), tolerance);
	}
}

TEST_F(Solver2DSpectralTest, DISABLED_periodic_tris_spectral_and_base_comparison)
{

	Probes probes;

	maxwell::Solver solver{
		buildModel(14, 1, Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(2.0)
			.setOrder(3)};

	maxwell::Solver solverSpectral{
		buildModel(14, 1, Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC, BdrCond::PEC, BdrCond::PMC, BdrCond::PEC),
		probes,
		buildGaussianInitialField(E, 0.1, fieldCenter, unitVec(Z)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(2.0)
			.setOrder(3)
			.setSpectralEO()};

	Vector zeroVec{solver.getField(E, Z).Size()};
	zeroVec = 0.0;
	double tolerance{1e-5};
	for (int i = 0; i < zeroVec.Size(); ++i)
	{
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E, Z).Elem(i) - solverSpectral.getField(E, Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, Y).Elem(i) - solverSpectral.getField(H, Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, X).Elem(i) - solverSpectral.getField(H, X).Elem(i), tolerance);
	}

	solver.run();
	solverSpectral.run();

	EXPECT_NEAR(solver.getFields().getNorml2(), solverSpectral.getFields().getNorml2(), 1e-4);
	for (int i = 0; i < zeroVec.Size(); ++i)
	{
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(E, Z).Elem(i) - solverSpectral.getField(E, Z).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, Y).Elem(i) - solverSpectral.getField(H, Y).Elem(i), tolerance);
		EXPECT_NEAR(zeroVec.Elem(i), solver.getField(H, X).Elem(i) - solverSpectral.getField(H, X).Elem(i), tolerance);
	}
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

TEST_F(Solver2DTest, periodic_tris)
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

TEST_F(Solver2DTest, periodic_quads)
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

TEST_F(Solver2DTest, DISABLED_pec_totalfieldin_1dot5D)
{
	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "4x4_Quadrilateral_1dot5D_IntBdr.mesh").c_str(), 1, 0)};
	mesh.UniformRefinement();
	// mesh.UniformRefinement();
	// mesh.UniformRefinement();
	GeomTagToBoundary attToBdr{{1, BdrCond::PEC}, {2, BdrCond::PMC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(5)};
	// probes.pointProbes = {
	// PointProbe{ E, Z, {0.5001, 0.5} },
	// PointProbe{ E, Z, {0.5, 0.5} },
	// PointProbe{ H, Y, {3.5, 0.5} },
	// PointProbe{ H, X, {3.5, 0.5} }
	// };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.5, 2.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			//.setCentered()
			.setFinalTime(6.0)
			.setOrder(3)};

	solver.run();

	//{
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().first, 1e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(3.5, solver.getPointProbe(2).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, 1e-3);
	//	EXPECT_NEAR(3.5, solver.getPointProbe(3).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, 1e-3);
	//}
}

TEST_F(Solver2DTest, DISABLED_pec_planewave)
{
	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "4x4_Quadrilateral_InnerSquare_IntBdr.mesh").c_str(), 1, 0)};
	// mesh.UniformRefinement();
	// mesh.UniformRefinement();
	// mesh.UniformRefinement();
	GeomTagToBoundary attToBdr{{1, BdrCond::PEC}, {2, BdrCond::PMC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(5)};
	// probes.pointProbes = {
	// PointProbe{ E, Z, {0.5001, 0.5} },
	// PointProbe{ E, Z, {0.5, 0.5} },
	// PointProbe{ H, Y, {3.5, 0.5} },
	// PointProbe{ H, X, {3.5, 0.5} }
	// };

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			//.setCentered()
			.setFinalTime(8.0)
			.setOrder(3)};

	solver.run();

	//{
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().first, 1e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(0).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, 1e-3);
	//}

	//{
	//	EXPECT_NEAR(3.5, solver.getPointProbe(2).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, 1e-3);
	//	EXPECT_NEAR(3.5, solver.getPointProbe(3).findFrameWithMax().first, 2e-1);
	//	EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, 1e-3);
	//}
}

TEST_F(Solver2DTest, DISABLED_sma_totalfieldinout_1dot5D)
{
	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "4x4_Quadrilateral_1dot5D_IntBdr.mesh").c_str(), 1, 0)};
	GeomTagToBoundary attToBdr{{1, BdrCond::SMA}, {2, BdrCond::PMC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(20)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 1.5, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(4.0)
			.setOrder(2)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_pec_centered_totalfieldin_longline_1dot5D)
{
	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "One_Element_Tall_Long_Line_TF_centered.mesh").c_str(), 1, 0)};
	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(20)};
	probes.pointProbes = {
		PointProbe{E, Z, {2.0, 0.125}},
		PointProbe{H, Y, {2.0, 0.125}},
		PointProbe{E, Z, {1.5, 0.125}},
		PointProbe{H, Y, {1.5, 0.125}}};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.2, 0.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setCentered()
			.setFinalTime(6.0)
			.setOrder(3)};

	solver.run();

	{
		auto frame{solver.getPointProbe(0).getFieldMovie()};
		auto expected_t = 2.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}
	{
		auto frame{solver.getPointProbe(1).getFieldMovie()};
		auto expected_t = 2.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}
	{
		auto frame{solver.getPointProbe(2).getFieldMovie()};
		auto expected_t = 5.5;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}
	{
		auto frame{solver.getPointProbe(3).getFieldMovie()};
		auto expected_t = 5.5;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}
}

TEST_F(Solver2DTest, DISABLED_pec_upwind_beam_totalfieldscatteredfield_in)
{
	Mesh mesh{Mesh::LoadFromFile((gmshMeshesFolder() + "2D_TF_Beam.msh").c_str(), 1, 0)};
	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(20)};
	probes.pointProbes = {
		PointProbe{E, Z, {0.0, 0.5}},
		PointProbe{E, Z, {1.0001, 0.5}},
		PointProbe{E, Z, {4.0, 0.5}},
		PointProbe{H, Y, {0.0, 0.5}},
		PointProbe{H, Y, {1.0001, 0.5}},
		PointProbe{H, Y, {4.0, 0.5}}};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.5, 1.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(2e-2)
			.setFinalTime(9.0)
			.setOrder(3)};

	solver.run();

	{
		auto frame{solver.getPointProbe(0).getFieldMovie()};
		auto expected_t = 1.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(0.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(0).getFieldMovie()};
		auto expected_t = 9.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(0.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(1).getFieldMovie()};
		auto expected_t = 2.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(2).getFieldMovie()};
		auto expected_t = 5.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(0.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(3).getFieldMovie()};
		auto expected_t = 9.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-2.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(4).getFieldMovie()};
		auto expected_t = 2.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(5).getFieldMovie()};
		auto expected_t = 5.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-2.0, f, 1e-2);
			}
		}
	}
}

TEST_F(Solver2DTest, DISABLED_pec_upwind_beam_totalfieldscatteredfield_inout)
{
	Mesh mesh{Mesh::LoadFromFile((gmshMeshesFolder() + "2D_TFSF_Beam.msh").c_str(), 1, 0)};
	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(20)};
	probes.pointProbes = {
		PointProbe{E, Z, {0.5, 0.5}},
		PointProbe{E, Z, {2.0, 0.5}},
		PointProbe{E, Z, {3.5, 0.5}},
		PointProbe{H, Y, {0.5, 0.5}},
		PointProbe{H, Y, {2.0, 0.5}},
		PointProbe{H, Y, {3.5, 0.5}}};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.5, 1.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(2e-2)
			.setFinalTime(9.0)
			.setOrder(3)};

	solver.run();

	{
		auto frame{solver.getPointProbe(0).findFrameWithMin()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}
	{
		auto frame{solver.getPointProbe(0).findFrameWithMax()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}

	{
		auto frame{solver.getPointProbe(1).getFieldMovie()};
		auto expected_t = 3.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(2).findFrameWithMin()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}
	{
		auto frame{solver.getPointProbe(2).findFrameWithMax()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}

	{
		auto frame{solver.getPointProbe(3).findFrameWithMin()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}
	{
		auto frame{solver.getPointProbe(3).findFrameWithMax()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}

	{
		auto frame{solver.getPointProbe(4).getFieldMovie()};
		auto expected_t = 3.0;
		for (const auto &[t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2)
			{
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{solver.getPointProbe(5).findFrameWithMin()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}
	{
		auto frame{solver.getPointProbe(5).findFrameWithMax()};
		EXPECT_NEAR(frame.second, 0.0, 1e-3);
	}
}

TEST_F(Solver2DTest, DISABLED_upwind_beam_totalfieldscatteredfield_in_fullface)
{
	Mesh mesh{Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "2D_DualSurface_FullFace_Beam.msh").c_str(), 1, 0)};
	// mesh.UniformRefinement();
	GeomTagToBoundary attToBdr{{2, BdrCond::PMC}, {3, BdrCond::PEC}};
	GeomTagToInteriorBoundary attToIntCond{{4, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, attToIntCond) };

	auto probes{buildProbesWithAnExportProbe(100)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.5, 2.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(9.0)
			.setOrder(3)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_upwind_beam_totalfieldscatteredfield_in_fullface_RtL)
{
	Mesh mesh{Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "2D_DualSurface_FullFace_Beam_RtL.msh").c_str(), 1, 0)};
	// mesh.UniformRefinement();
	GeomTagToBoundary attToBdr{{2, BdrCond::PMC}, {3, BdrCond::PEC}};
	GeomTagToInteriorBoundary attToIntCond{{4, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, attToIntCond) };

	auto probes{buildProbesWithAnExportProbe(100)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.5, 4.0, unitVec(Z), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(9.0)
			.setOrder(3)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_pec_upwind_totalfieldin_square_1dot5D)
{
	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "4x4_Quadrilateral_InnerSquare_IntBdr.mesh").c_str(), 1, 0)};
	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(25)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 1.5, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(10.0)
			.setOrder(3)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_pec_upwind_totalfieldin_square_1dot5D_rotated45)
{
	Mesh mesh{Mesh::LoadFromFile((mfemMeshes2DFolder() + "4x4_Quadrilateral_InnerSquare_IntBdr.mesh").c_str(), 1, 0)};
	GeomTagToBoundary attToBdr{{1, BdrCond::PMC}, {2, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(25)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 1.0, unitVec(Z), Vector{{1.0 / sqrt(2.0), 1.0 / sqrt(2.0), 0.0}}),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(7.0)
			.setOrder(3)};

	solver.run();
}

// TEST_F(Solver2DTest, DISABLED_quadraticMesh)
//{
//	Mesh mesh = Mesh::LoadFromFile("./testData/star-q2.mesh", 1, 0);
//	auto fec = std::make_unique<DG_FECollection>(4, 2, BasisType::GaussLobatto);
//	auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get());
//
//	Model model = Model(mesh, GeomTagToMaterial{}, GeomTagToBoundary{});
//
//	maxwell::Solver solver{
//		model,
//		buildProbesWithAnExportProbe(),
//		buildGaussianInitialField(E, Z, 0.4, mfem::Vector({-0.02566,0.03028})),
//		SolverOptions{}
//		.setTimeStep(5e-4)
//		.setFinalTime(2.0)
//		.setOrder(3)
//	};
//
//	auto normOld{ solver.getFields().getNorml2() };
//	solver.run();
//
//	double tolerance{ 1e-2 };
//	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);
//
// }

// TEST_F(Solver2DTest, DISABLED_periodic_strong) //TODO ADD ENERGY CHECK
//{
//	Mesh mesh2D = Mesh::MakeCartesian2D(21, 3, Element::Type::QUADRILATERAL);
//	Vector periodic({ 0.0, 1.0 });
//	std::vector<Vector> trans;
//	trans.push_back(periodic);
//	Mesh mesh2DPer = Mesh::MakePeriodic(mesh2D, mesh2D.CreatePeriodicVertexMapping(trans));
//
//	maxwell::Solver::Options opts;
//	opts.evolution = FiniteElementEvolution::Options();
//	opts.evolution.disForm = DisForm::Strong;
//
//	Model model = Model(mesh2DPer, GeomTagToMaterial(), GeomTagToBoundary());
//
//	Probes probes;
//	probes.addExporterProbeToCollection(ExporterProbe());
//	probes.vis_steps = 20;
//
//	Sources sources;
//	sources.addSourceToVector(Source(model, E, X, 1.0, 10.0, Vector({ 0.0, 0.0 })));
//
//	maxwell::Solver solver(
//		model,
//		probes,
//		sources,
//		opts);
//
//	solver.run();
//
// }
// TEST_F(Solver2DTest, DISABLED_centered_NC_MESH) //TODO ADD ENERGY CHECK
//{
//	/*The purpose of this test is to verify the functionality of the Maxwell Solver when using
//	a centered type flux. A non-conforming mesh is loaded to test MFEM functionalities on the code.
//
//	First, all required parts for constructing a solver are declared, Model, Sources, Probes and Options.
//	A single 2D Gaussian on X and Y is declared along Ez.
//
//	Then, the Solver object is constructed using said parts, with its mesh being two-dimensional mixed
//	with triangles and squares.
//	The field along Ez is extracted before and after the solver calls its run() method and evolves the
//	problem. This test verifies that, for this mesh, after two seconds and nine hundred twenty
//	miliseconds, the problem reaches a new peak in field Ez and the maximum value in Ez is not
//	higher than the initial value.*/
//
//	const char* mesh_file = "star-mixed.mesh";
//	Mesh mesh(mesh_file);
//	mesh.UniformRefinement();
//	Model model = Model(mesh, GeomTagToMaterial(), GeomTagToBoundary());
//
//	Probes probes;
//	//probes.addExporterProbeToCollection(ExporterProbe());
//	//probes.vis_steps = 20;
//
//	Sources sources;
//	sources.addSourceToVector(Source(model, E, Z, 2.0, 20.0, Vector({ 0.0, 0.0 })));
//
//	maxwell::Solver::Options solverOpts = buildDefaultSolverOpts(2.92);
//	solverOpts.evolution.fluxType = FluxType::Centered;
//
//	maxwell::Solver solver(model, probes, sources, solverOpts);
//
//	GridFunction eOld = solver.getFieldInDirection(E, Z);
//	solver.run();
//	GridFunction eNew = solver.getFieldInDirection(E, Z);
//
//	EXPECT_GT(eOld.Max(), eNew.Max());
// }


TEST_F(Solver2DTest, interiorPEC_sma_boundaries)
{
	Mesh mesh{Mesh::LoadFromFile((gmshMeshesFolder() + "InteriorPEC2D.msh").c_str(), 1, 0)};
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

TEST_F(Solver2DTest, interiorBoundary_TotalFieldIn)
{
	auto mesh{
		Mesh::LoadFromFile(
			(mfemMeshes2DFolder() + "intbdr_two_quads.mesh").c_str(), 1, 0)};
	// mesh.UniformRefinement();
	GeomTagToBoundary attToBdr{{1, BdrCond::PEC}, {2, BdrCond::PMC}};

	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };
	auto probes{buildProbesWithAnExportProbe(80)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.3, 0.5, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setCentered()
			.setFinalTime(2.0)
			.setOrder(2)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_box_with_Gmsh)
{
	auto mesh = Mesh::LoadFromFile((gmshMeshesFolder() + "test.msh").c_str(), 1, 0);
	auto fec = std::make_unique<DG_FECollection>(1, 2, BasisType::GaussLobatto);
	auto fes = std::make_unique<FiniteElementSpace>(&mesh, fec.get(), 1, 0);
	auto model = Model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}));

	auto probes{buildProbesWithAnExportProbe(100)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianInitialField(E, 0.1, mfem::Vector({0.5, 0.5})),
		SolverOptions{}
			.setTimeStep(5e-4)
			.setFinalTime(2.0)
			.setOrder(4)};
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

TEST_F(Solver2DTest, DISABLED_upwind_box_totalfieldscatteredfield_inout_circleinside)
{
	Mesh mesh{Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "2D_TFSF_Box_Circle.msh").c_str(), 1, 0, true)};
	GeomTagToBoundary attToBdr{{2, BdrCond::SMA}, {3, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(20)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 1.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(5.0)
			.setOrder(2)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_upwind_beam_totalfieldscatteredfield_in_intbdr_fss)
{
	Mesh mesh{Mesh::LoadFromFile((gmshMeshesFolder() + "2D_TF_FSS.msh").c_str(), 1, 0, true)};
	GeomTagToBoundary attToBdr{{2, BdrCond::SMA}, {3, BdrCond::PMC}};
	GeomTagToInteriorBoundary attToIntCond{{4, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, attToIntCond) };

	auto probes{buildProbesWithAnExportProbe(10)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 2.5, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(10.0)
			.setOrder(3)};

	solver.run();
}

TEST_F(Solver2DTest, DISABLED_upwind_box_totalfieldscatteredfield_inout_circle_w_circles)
{
	Mesh mesh{Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "2D_TFSF_Circle_of_circles.msh").c_str(), 1, 0, true)};
	GeomTagToBoundary attToBdr{{2, BdrCond::SMA}, {3, BdrCond::PEC}};
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, GeomTagToInteriorBoundary{}) };

	auto probes{buildProbesWithAnExportProbe(20)};

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.5, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-3)
			.setFinalTime(30.0)
			.setOrder(3)};

	solver.run();
}