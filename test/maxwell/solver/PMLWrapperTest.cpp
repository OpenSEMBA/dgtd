#include <gtest/gtest.h>

#include "evolution/Fields.h"
#include "solver/SolverExtension.h"

using namespace maxwell;
using namespace mfem;

namespace {

Mesh make1DPMLShellMesh(int n = 5)
{
	auto m = Mesh::MakeCartesian1D(n);
	m.GetElement(0)->SetAttribute(2);
	m.GetElement(n - 1)->SetAttribute(2);
	m.SetAttributes();
	return m;
}

}

class PMLWrapperTest : public ::testing::Test {
};

TEST_F(PMLWrapperTest, buildsFromTaggedVolumes)
{
	auto mesh = make1DPMLShellMesh();
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(1.0, 1.0, 0.0));
	Model model(mesh, matInfo);
	GeomTagToPMLRegion pml_regions;
	pml_regions.emplace(2, PMLRegionProperties{});
	model.setPMLRegions(pml_regions);

	mfem::DG_FECollection fec(1, mesh.Dimension());
	mfem::ParFiniteElementSpace parent_fes(&model.getMesh(), &fec);

	PMLWrapper wrapper(model, parent_fes);
	PMLRegionState state;
	state.init(wrapper.getStateSize());

	ASSERT_NE(nullptr, wrapper.getConstSubMesh());
	EXPECT_EQ(wrapper.getStateSize(), state.auxiliary_state.Size());
	EXPECT_NEAR(0.2, wrapper.getGeometry().thickness_minus[0], 1e-12);
	EXPECT_NEAR(0.2, wrapper.getGeometry().thickness_plus[0], 1e-12);
	EXPECT_GT(wrapper.getSubToParentVDofMap().Size(), 0);

	const int local_parent_elements = wrapper.getConstSubMesh()->GetParentElementIDMap().Size();
	int global_parent_elements = 0;
	MPI_Allreduce(&local_parent_elements, &global_parent_elements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	EXPECT_EQ(2, global_parent_elements);

	for (int i = 0; i < wrapper.getSubToParentVDofMap().Size(); ++i) {
		EXPECT_LT(wrapper.getSubToParentVDofMap()[i], parent_fes.GetVSize());
		EXPECT_LE(0, wrapper.getSubToParentVDofMap()[i]);
	}
}

TEST_F(PMLWrapperTest, implicitCorrectionDampsOnlyPMLDofs)
{
	auto mesh = make1DPMLShellMesh(9);
	GeomTagToMaterialInfo matInfo;
	matInfo.gt2m.emplace(1, Material(1.0, 1.0, 0.0));
	matInfo.gt2m.emplace(2, Material(1.0, 1.0, 0.0));
	Model model(mesh, matInfo);
	GeomTagToPMLRegion pml_regions;
	pml_regions.emplace(2, PMLRegionProperties{});
	model.setPMLRegions(pml_regions);

	mfem::DG_FECollection fec(1, mesh.Dimension());
	mfem::ParFiniteElementSpace parent_fes(&model.getMesh(), &fec);
	Fields<mfem::ParFiniteElementSpace, mfem::ParGridFunction> parent_fields(parent_fes);
	for (int dof = 0; dof < parent_fes.GetNDofs(); ++dof) {
		parent_fields.get(E, Y)[dof] = 1.0;
		parent_fields.get(H, Z)[dof] = 1.0;
	}

	PMLWrapper wrapper(model, parent_fes);
	PMLRegionState state;
	state.init(wrapper.getStateSize());
	wrapper.applyImplicitCorrection(0.1, parent_fields, state);

	std::vector<bool> is_pml_vdof(parent_fes.GetNDofs(), false);
	for (int i = 0; i < wrapper.getSubToParentVDofMap().Size(); ++i) {
		const int mapped = wrapper.getSubToParentVDofMap()[i] >= 0
			? wrapper.getSubToParentVDofMap()[i]
			: -1 - wrapper.getSubToParentVDofMap()[i];
		is_pml_vdof[mapped] = true;
	}

	bool found_damped = false;
	bool found_interior = false;
	for (int dof = 0; dof < parent_fes.GetNDofs(); ++dof) {
		if (is_pml_vdof[dof]) {
			EXPECT_LT(parent_fields.get(E, Y)[dof], 1.0);
			EXPECT_LT(parent_fields.get(H, Z)[dof], 1.0);
			found_damped = true;
		} else {
			EXPECT_DOUBLE_EQ(1.0, parent_fields.get(E, Y)[dof]);
			EXPECT_DOUBLE_EQ(1.0, parent_fields.get(H, Z)[dof]);
			found_interior = true;
		}
	}

	EXPECT_TRUE(found_damped);
	EXPECT_TRUE(found_interior);
	EXPECT_GT(wrapper.getSigmaMax(), 0.0);
	EXPECT_GT(wrapper.getMinNormalElementSize(), 0.0);
}