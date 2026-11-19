#include "components/Model.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include <unistd.h>

namespace {

using namespace maxwell;

class LorentzGridFunctionTest : public ::testing::Test {
protected:
	std::filesystem::path dir;

	void SetUp() override
	{
		dir = std::filesystem::temp_directory_path() /
			("dgtd_lorentz_gf_" + std::to_string(static_cast<long long>(::getpid())));
		std::filesystem::create_directories(dir);
	}

	void TearDown() override
	{
		std::error_code ec;
		std::filesystem::remove_all(dir, ec);
	}
};

} // namespace

TEST_F(LorentzGridFunctionTest, LoadsAndDistributesToParallelSpace)
{
	mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
		2, 2, 2, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
	mfem::DG_FECollection fec(2, 3, mfem::BasisType::GaussLobatto);
	mfem::FiniteElementSpace serial_fes(&mesh, &fec);
	mfem::GridFunction serial_gf(&serial_fes);
	mfem::FunctionCoefficient coefficient(
		[](const mfem::Vector& x) { return 2.0 + x(0) + 0.5 * x(1); });
	serial_gf.ProjectCoefficient(coefficient, mfem::ProjectType::ELEMENT_L2);

	const std::filesystem::path gf_path = dir / "omega_p.gf";
	serial_gf.Save(gf_path.c_str(), 12);

	GeomTagToMaterialInfo materials;
	materials.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	GeomTagToBoundaryInfo boundaries;
	Model model(mesh, materials, boundaries);

	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = 1.0;
	pole.omega_p_file = gf_path.string();
	model.setLorentzProperties({pole});

	mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
	mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
	model.loadLorentzGridFunctions(pfes);

	const LorentzProperties* loaded = model.findLorentz(1);
	ASSERT_NE(loaded, nullptr);
	ASSERT_NE(loaded->omega_p_gf, nullptr);
	EXPECT_GT(loaded->omega_p, 0.0);

	mfem::Array<int> serial_dofs, local_dofs;
	for (int e = 0; e < pmesh.GetNE(); ++e) {
		const int global = static_cast<int>(pmesh.GetGlobalElementNum(e));
		serial_fes.GetElementDofs(global, serial_dofs);
		pfes.GetElementDofs(e, local_dofs);
		ASSERT_EQ(serial_dofs.Size(), local_dofs.Size());
		for (int i = 0; i < local_dofs.Size(); ++i) {
			EXPECT_NEAR(
				(*loaded->omega_p_gf)(local_dofs[i]),
				serial_gf(serial_dofs[i]), 1.0e-12);
		}
	}
}

TEST_F(LorentzGridFunctionTest, RejectsAGridFunctionOnADifferentMesh)
{
	mfem::Mesh coarse = mfem::Mesh::MakeCartesian3D(
		1, 1, 1, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
	mfem::DG_FECollection fec(2, 3, mfem::BasisType::GaussLobatto);
	mfem::FiniteElementSpace serial_fes(&coarse, &fec);
	mfem::GridFunction serial_gf(&serial_fes);
	serial_gf = 1.0;
	const std::filesystem::path gf_path = dir / "wrong_mesh.gf";
	serial_gf.Save(gf_path.c_str(), 12);

	mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
		2, 2, 2, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
	GeomTagToMaterialInfo materials;
	materials.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	GeomTagToBoundaryInfo boundaries;
	Model model(mesh, materials, boundaries);
	LorentzProperties pole;
	pole.geom_tag = 1;
	pole.eps_inf = 1.0;
	pole.omega_p_file = gf_path.string();
	model.setLorentzProperties({pole});

	mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
	mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
	EXPECT_THROW(model.loadLorentzGridFunctions(pfes), std::runtime_error);
}
