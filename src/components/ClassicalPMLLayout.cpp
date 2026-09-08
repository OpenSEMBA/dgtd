#include "ClassicalPMLLayout.h"

namespace maxwell {

ClassicalPMLLayout::ClassicalPMLLayout(
	int ndofs, const std::vector<PMLProperties>& regions, int mesh_dim)
	: ndofs_(ndofs)
{
	if (regions.empty() || ndofs_ <= 0) {
		return;
	}
	if (mesh_dim < 1 || mesh_dim > 3) {
		throw std::runtime_error("ClassicalPMLLayout: invalid mesh dimension.");
	}

	// Multi-axis blocks (e.g. active_axes: ["X","Y"]) are allowed: biaxial ADE is
	// the superposition of per-axis uniaxial stacks on the same elements. Depth
	// ρ_s / L_s remain axis-aligned (planar vacuum–PML interfaces); the user must
	// ensure the mesh geometry matches that model (centered box/ring, etc.).
	for (const auto& props : regions) {
		for (Direction d : props.active_axes) {
			if (d < 0 || d >= mesh_dim) {
				throw std::runtime_error(
					"ClassicalPMLLayout: active_axes exceeds mesh dimension.");
			}
			if (dir_to_slot_[d] < 0) {
				dir_to_slot_[d] = static_cast<int>(stretch_dirs_.size());
				stretch_dirs_.push_back(d);
			}
		}
	}

	n_aux_ = 2 * static_cast<int>(stretch_dirs_.size()) * ndofs_;
}

bool ClassicalPMLLayout::isActiveDirection(Direction d) const
{
	return d >= X && d <= Z && dir_to_slot_[d] >= 0;
}

int ClassicalPMLLayout::stretchSlot(Direction d) const
{
	if (!isActiveDirection(d)) {
		throw std::runtime_error("ClassicalPMLLayout: direction is not active.");
	}
	return dir_to_slot_[d];
}

int ClassicalPMLLayout::jOffset(Direction stretch_d) const
{
	return 6 * ndofs_ + 2 * stretchSlot(stretch_d) * ndofs_;
}

int ClassicalPMLLayout::mOffset(Direction stretch_d) const
{
	return 6 * ndofs_ + (2 * stretchSlot(stretch_d) + 1) * ndofs_;
}

} // namespace maxwell
