#include "PMLAuxLayout.h"

#include "Model.h"

#include <algorithm>
#include <set>

namespace maxwell {

namespace {

std::set<Direction> collectStretchDirections(
	const std::vector<PMLProperties>& pml_props, int mesh_dim)
{
	std::set<Direction> dirs;
	for (const auto& props : pml_props) {
		for (Direction d : props.active_axes) {
			if (d >= 0 && d < mesh_dim) {
				dirs.insert(d);
			}
		}
	}
	return dirs;
}

} // namespace

PMLAuxLayout::PMLAuxLayout(
	const mfem::ParFiniteElementSpace& fes, const std::vector<PMLProperties>& pml_props)
{
	ndofs_ = fes.GetNDofs();
	field_block_size_ = 6 * ndofs_;
	if (pml_props.empty()) {
		return;
	}

	const int mesh_dim = fes.GetMesh()->Dimension();
	const std::set<Direction> dirs = collectStretchDirections(pml_props, mesh_dim);
	stretch_dirs_.assign(dirs.begin(), dirs.end());
	dir_to_slot_.assign(3, -1);
	for (int i = 0; i < static_cast<int>(stretch_dirs_.size()); ++i) {
		dir_to_slot_[stretch_dirs_[i]] = i;
	}
	n_aux_ = auxBlockSize() * static_cast<int>(stretch_dirs_.size());
}

int PMLAuxLayout::psiEOffset(Direction stretch_d, Direction h_comp) const
{
	const int slot = stretchSlot(stretch_d);
	if (slot < 0 || h_comp < X || h_comp > Z) {
		return -1;
	}
	return field_block_size_ + slot * auxBlockSize() + static_cast<int>(h_comp) * ndofs_;
}

int PMLAuxLayout::psiHOffset(Direction stretch_d, Direction e_comp) const
{
	const int slot = stretchSlot(stretch_d);
	if (slot < 0 || e_comp < X || e_comp > Z) {
		return -1;
	}
	return field_block_size_ + slot * auxBlockSize() + 3 * ndofs_ +
	       static_cast<int>(e_comp) * ndofs_;
}

bool PMLAuxLayout::isActiveDirection(Direction d) const
{
	return d >= 0 && d < 3 && dir_to_slot_[d] >= 0;
}

int PMLAuxLayout::stretchSlot(Direction d) const
{
	if (d < 0 || d >= 3) {
		return -1;
	}
	return dir_to_slot_[d];
}

int computePMLAuxSize(const Model& model, int ndofs, int mesh_dim)
{
	if (!model.hasPML()) {
		return 0;
	}
	const std::set<Direction> dirs =
		collectStretchDirections(model.getPMLProperties(), mesh_dim);
	return 6 * ndofs * static_cast<int>(dirs.size());
}

} // namespace maxwell
