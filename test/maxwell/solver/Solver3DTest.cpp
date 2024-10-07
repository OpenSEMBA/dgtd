#include <gtest/gtest.h>

#include "ProbeFixtures.h"
#include "SourceFixtures.h"

#include "solver/Solver.h"

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

TEST_F(Solver3DTest, pec_1dot5D)
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

TEST_F(Solver3DTest, periodic_cube_hexa)
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
				Gaussian{0.1}, 
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
	maxwell::Solver solver{
	buildModel(
		12,    4,   4, Element::Type::HEXAHEDRON,
		30.0, 15.0, 12.0,
		BdrCond::PEC,BdrCond::PMC,BdrCond::SMA,
		BdrCond::PMC,BdrCond::SMA,BdrCond::PEC),
	buildProbesEmpty(),
	buildGaussianInitialField(
		E, 3.0,
		Source::Position({15.0,7.5,6.0}),
		unitVec(Z)
	),
	SolverOptions{}
		.setTimeStep(5e-1)
		.setFinalTime(5.0)
		.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DTest, interiorPEC_sma_boundaries)
{
	Mesh mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "InteriorPEC3D.msh").c_str(),1,0)};
	GeomTagToBoundary attToBdr{ {2,BdrCond::PEC},{3,BdrCond::PMC}, {4,BdrCond::SMA } };
	GeomTagToInteriorBoundary attToIntBdr{ {5,BdrCond::PEC} };
	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(attToBdr, attToIntBdr) };

	auto probes{ buildProbesWithAnExportProbe() };

	maxwell::Solver solver{
		model,
		probes,
		buildPlanewaveInitialField(
			Gaussian{7.5},
			Source::Position({ 0.0 }), // center
			Source::Polarization(unitVec(Z)), // e polarization
			mfem::Vector({1.0, 0.0, 0.0}) // propagation direction
	),
		SolverOptions{}
			.setTimeStep(1.0)
			.setFinalTime(30.0)
			.setOrder(2)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

}

TEST_F(Solver3DTest, DISABLED_minimal_tetra)
{
	auto probes{ buildProbesWithAnExportProbe() };
	Mesh mesh{ Mesh::LoadFromFile((gmshMeshesFolder() + "twotetras.msh").c_str(),1,0)};

	GeomTagToBoundary attToBdr{ {1, BdrCond::PEC},{2, BdrCond::PMC},{3, BdrCond::PMC},{4, BdrCond::PMC},{5, BdrCond::PMC},{6, BdrCond::PEC} };

	Model model{ mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(GeomTagToBoundary{}, GeomTagToInteriorBoundary{}) };

	maxwell::Solver solver{
	model,
	probes,
	buildPlanewaveInitialField(
		Gaussian{0.16},
		Source::Position({ 0.7, 0.0, 0.0 }), // center
		Source::Polarization(unitVec(Z)), // e polarization
		mfem::Vector(unitVec(X)) // propagation direction
	),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setFinalTime(2.0)
		.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, sma_upwind_hex_1dot5D)
{

	Probes probes{ buildProbesWithAnExportProbe(200) };
	//probes.pointProbes = {
	//	PointProbe{E, Z, {0.0, 0.5, 0.5}},
	//	PointProbe{E, Z, {1.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {0.0, 0.5, 0.5}},
	//	PointProbe{H, Y, {1.0, 0.5, 0.5}}
	//};

	maxwell::Solver solver{
		buildModel(
			10, 2, 2, Element::Type::HEXAHEDRON,
			1.0, 0.4, 0.4,
			BdrCond::PEC, BdrCond::PMC, BdrCond::SMA,
			BdrCond::PMC, BdrCond::SMA, BdrCond::PEC),
			probes,
			buildGaussianInitialField(
				E, 0.1,
				Source::Position({ 0.5,0.5,0.5 }),
				unitVec(Z)
			),
			SolverOptions{}
			.setTimeStep(1e-4)
			.setFinalTime(2.0)
			.setOrder(3)
	};

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//EXPECT_NEAR(0.0, solver.getPointProbe(0).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(0.0, solver.getPointProbe(1).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(2).findFrameWithMax().second, tolerance);
	//EXPECT_NEAR(1.0, solver.getPointProbe(3).findFrameWithMax().second, tolerance);
}

TEST_F(Solver3DTest, DISABLED_pec_centered_hexa_totalfieldin)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {3.0, 0.5, 0.5}},
		PointProbe{E, Z, {5.0, 0.5, 0.5}},
		PointProbe{H, Y, {3.0, 0.5, 0.5}},
		PointProbe{H, Y, {5.0, 0.5, 0.5}}
	};
	auto mesh{ Mesh::LoadFromFile((mfemMeshes3DFolder() + "beam_hex_totalfieldin.mesh").c_str(), 1, 0) };
	GeomTagToBoundary att2bdr{ {1, BdrCond::PMC}, {2, BdrCond::PEC}, {3, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, GeomTagToInteriorBoundary()));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(8.0)
			.setOrder(3)
	};

	solver.run();
	
	{
		auto frame{ solver.getPointProbe(0).getFieldMovie()};
		auto expected_t = 6.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(1).getFieldMovie() };
		auto expected_t = 8.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(2).getFieldMovie() };
		auto expected_t = 6.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}

	{
		auto frame{ solver.getPointProbe(3).getFieldMovie() };
		auto expected_t = 8.0;
		for (const auto& [t, f] : frame)
		{
			if (abs(expected_t - t) <= 1e-2) {
				EXPECT_NEAR(-1.0, f, 1e-2);
			}
		}
	}

}

TEST_F(Solver3DTest, DISABLED_pec_upwind_box_totalfieldscatteredfield)
{
	auto probes{ buildProbesWithAnExportProbe(30) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_MinimalistBox.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::SMA} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, GeomTagToInteriorBoundary()));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.5, 5.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(5e-3)
			.setFinalTime(10.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_pec_centered_beam_totalfieldscatteredfield)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {2.0, 0.5, 0.5}},
		PointProbe{E, Z, {5.0, 0.5, 0.5}},
		PointProbe{H, Y, {2.0, 0.5, 0.5}},
		PointProbe{H, Y, {5.0, 0.5, 0.5}}
	};
	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PMC}, {1, BdrCond::PEC}, {3, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, GeomTagToInteriorBoundary()));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_upwind_beam_totalfieldscatteredfield_inout)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	probes.pointProbes = {
		PointProbe{E, Z, {2.0, 0.5, 0.5}},
		PointProbe{E, Z, {5.0, 0.5, 0.5}},
		PointProbe{H, Y, {2.0, 0.5, 0.5}},
		PointProbe{H, Y, {5.0, 0.5, 0.5}}
	};
	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PMC}, {1, BdrCond::PEC}, {3, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, GeomTagToInteriorBoundary()));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_dualintbdr_upwind_beam_totalfieldscatteredfield_in)
{
	auto probes{ buildProbesWithAnExportProbe(10) };
	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_DualSurface_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC}, {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, GeomTagToInteriorBoundary()));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(0.7, 2.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(10.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_pec_centered_innerbox_totalfieldinout)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TFSF_Box.msh").c_str(), 1, 0, true) };
	mesh.UniformRefinement();
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, GeomTagToInteriorBoundary()));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(10.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_centered_beam_totalfieldscatteredfield_inout_intbdr)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorBoundary att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, att2IntCond));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_centered_beam_totalfieldscatteredfield_inout_intbdr_RtL)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam_RtL.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorBoundary att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, att2IntCond));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 8.0, unitVec(Z), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setCentered()
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_upwind_beam_totalfieldscatteredfield_inout_intbdr)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorBoundary att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, att2IntCond));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 3.0, unitVec(Z), unitVec(X)),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}

TEST_F(Solver3DTest, DISABLED_upwind_beam_totalfieldscatteredfield_inout_intbdr_RtL)
{
	auto probes{ buildProbesWithAnExportProbe(10) };

	auto mesh{ Mesh::LoadFromFileNoBdrFix((gmshMeshesFolder() + "3D_TF_IntBdr_Beam_RtL.msh").c_str(), 1, 0, true) };
	GeomTagToBoundary att2bdr{ {2, BdrCond::PEC}, {3, BdrCond::PMC}, {4, BdrCond::PEC} };
	GeomTagToInteriorBoundary att2IntCond{ {5, BdrCond::PEC} };
	Model model(mesh, GeomTagToMaterialInfo(), GeomTagToBoundaryInfo(att2bdr, att2IntCond));

	maxwell::Solver solver{
		model,
		probes,
		buildGaussianPlanewave(1.0, 8.0, unitVec(Z), Vector{{-1.0, 0.0, 0.0}}),
		SolverOptions{}
			.setTimeStep(1e-2)
			.setFinalTime(20.0)
			.setOrder(3)
	};

	solver.run();

}