#include "math/Function.h"
#include "math/GridFunctionProfile.h"
#include "math/PhysicalConstants.h"
#include "components/LorentzProperties.h"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include <unistd.h>

namespace {

using namespace maxwell;

struct GridFunctionFixture {
	std::filesystem::path dir;
	std::filesystem::path mesh_path;
	std::filesystem::path gf_path;
	double scale = 0.0;
	double reference_value = 0.0;

	GridFunctionFixture()
	{
		dir = std::filesystem::temp_directory_path() /
			("dgtd_gridfunction_profile_" + std::to_string(::getpid()));
		std::filesystem::create_directories(dir);
		mesh_path = dir / "mesh.mesh";
		gf_path = dir / "omega_p.gf";

		mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(
			4, 2, mfem::Element::QUADRILATERAL, true, 1.0, 1.0);
		mesh.Save(mesh_path.string());

		mfem::H1_FECollection fec(1, 2);
		mfem::FiniteElementSpace fes(&mesh, &fec);
		mfem::GridFunction gf(&fes);
		for (int i = 0; i < fes.GetNDofs(); ++i) {
			const double* vertex = mesh.GetVertex(i);
			gf(i) = 1.0e9 * (1.0 + vertex[0] + 0.5 * vertex[1]);
		}
		gf.Save(gf_path.string().c_str());

		scale = 1.0 / physicalConstants::speedOfLight_SI;
		reference_value = 1.0e9 * (1.0 + 0.5 + 0.25); // (x,y) = (0.5,0.5)
	}

	~GridFunctionFixture()
	{
		std::error_code ec;
		std::filesystem::remove_all(dir, ec);
	}
};

} // namespace

TEST(GridFunctionProfile, EvaluatesSourceMeshAtInteriorPoints)
{
	GridFunctionFixture fx;
	GridFunctionProfile profile(
		fx.mesh_path.string(), fx.gf_path.string(), true, fx.scale, 0.0);

	EXPECT_EQ(profile.dimension(), 3); // axisymmetric reduction consumes 3D points
	EXPECT_TRUE(profile.axisymmetric());
	EXPECT_NEAR(profile.maxValue(), 1.0e9 * 2.5, 1.0);
	EXPECT_GT(profile.minValue(), 0.0);

	mfem::Vector p(3);
	p(0) = 0.5;
	p(1) = 0.5;
	p(2) = 0.0;
	EXPECT_NEAR(profile.eval(p), fx.reference_value * fx.scale, 1.0e-9);
}

TEST(GridFunctionProfile, AxisymmetricReductionUsesRadius)
{
	GridFunctionFixture fx;
	GridFunctionProfile profile(
		fx.mesh_path.string(), fx.gf_path.string(), true, fx.scale, 0.0);

	mfem::Vector on_plane(3);
	on_plane(0) = 0.5;
	on_plane(1) = 0.25;
	on_plane(2) = 0.0;

	mfem::Vector rotated(3);
	rotated(0) = 0.5;
	rotated(1) = 0.25 * std::cos(0.9);
	rotated(2) = 0.25 * std::sin(0.9);

	EXPECT_NEAR(profile.eval(on_plane), profile.eval(rotated), 1.0e-12);
}

TEST(GridFunctionProfile, OutsideValueIsReturnedOutsideTheSourceMesh)
{
	GridFunctionFixture fx;
	GridFunctionProfile profile(
		fx.mesh_path.string(), fx.gf_path.string(), true, fx.scale, 7.0);

	mfem::Vector outside(3);
	outside(0) = 100.0;
	outside(1) = 0.0;
	outside(2) = 0.0;
	EXPECT_NEAR(profile.eval(outside), 7.0 * fx.scale, 1.0e-15);
}

TEST(GridFunctionProfile, NearestFallbackCoversBoundaryMisses)
{
	GridFunctionFixture fx;
	GridFunctionProfile profile(
		fx.mesh_path.string(), fx.gf_path.string(), true, fx.scale, 0.0, true);

	// Just outside the source square, close to the (1, 0.1) boundary.
	mfem::Vector near_outside(3);
	near_outside(0) = 1.02;
	near_outside(1) = 0.1;
	near_outside(2) = 0.0;
	const double value = profile.eval(near_outside);
	EXPECT_GT(value, 0.0);
	// nearest element is [0.75,1]x[0,0.5]; its center is (0.875, 0.25)
	EXPECT_NEAR(value, 1.0e9 * (1.0 + 0.875 + 0.5 * 0.25) * fx.scale, 1.0e-9);
}

TEST(GridFunctionProfile, ParseLorentzObjectReadsGridFunctionRates)
{
	GridFunctionFixture fx;
	const std::string case_path = (fx.dir / "case.json").string();

	nlohmann::json lorentz{
		{"eps_inf", 1.0},
		{"f_p", {
			{"type", "gridfunction"},
			{"mesh", fx.mesh_path.string()},
			{"file", fx.gf_path.string()},
			{"axisymmetric", true}
		}},
		{"f_1", 0.0},
		{"gamma", 1.0e8}
	};

	auto pole = parseLorentzObject(lorentz, case_path);
	EXPECT_TRUE(pole.hasSpatialProfile());
	EXPECT_NE(pole.omega_p_profile, nullptr);
	EXPECT_EQ(pole.gamma_profile, nullptr);
	EXPECT_GT(pole.omega_p, 0.0);

	mfem::Vector p(3);
	p(0) = 0.5;
	p(1) = 0.5;
	p(2) = 0.0;
	EXPECT_NEAR(
		pole.evalOmegaP(p),
		fx.reference_value / physicalConstants::speedOfLight_SI,
		1.0e-9);
}

TEST(GridFunctionProfile, RejectsGridFunctionWithoutPositiveValues)
{
	GridFunctionFixture fx;
	mfem::Mesh mesh = mfem::Mesh::LoadFromFile(fx.mesh_path.string(), 1, 0, true);
	mfem::H1_FECollection fec(1, 2);
	mfem::FiniteElementSpace fes(&mesh, &fec);
	mfem::GridFunction gf(&fes);
	gf = 0.0;
	const std::filesystem::path zero_path = fx.dir / "zero.gf";
	gf.Save(zero_path.string().c_str());

	nlohmann::json lorentz{
		{"eps_inf", 1.0},
		{"f_p", {
			{"type", "gridfunction"},
			{"mesh", fx.mesh_path.string()},
			{"file", zero_path.string()}
		}},
		{"f_1", 0.0},
		{"gamma", 0.0}
	};
	EXPECT_THROW(parseLorentzObject(lorentz, ""), std::runtime_error);
}
