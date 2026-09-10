#include "SCPMLLayout.h"

namespace maxwell {

SCPMLLayout::SCPMLLayout(
	int ndofs, const std::vector<PMLProperties>& regions, int mesh_dim)
	: ndofs_(ndofs)
{
	if (regions.empty() || ndofs_ <= 0) {
		return;
	}
	if (mesh_dim < 1 || mesh_dim > 3) {
		throw std::runtime_error("SCPMLLayout: invalid mesh dimension.");
	}

	bool any_axis = false;
	for (const auto& props : regions) {
		for (Direction d : props.active_axes) {
			if (d < 0 || d >= mesh_dim) {
				throw std::runtime_error(
					"SCPMLLayout: active_axes exceeds mesh dimension.");
			}
			any_axis = true;
		}
	}
	if (any_axis) {
		n_aux_ = 6 * ndofs_;
	}
}

int SCPMLLayout::pEOffset(Direction comp) const
{
	if (comp < X || comp > Z) {
		throw std::runtime_error("SCPMLLayout: invalid P_E component.");
	}
	return 6 * ndofs_ + comp * ndofs_;
}

int SCPMLLayout::pHOffset(Direction comp) const
{
	if (comp < X || comp > Z) {
		throw std::runtime_error("SCPMLLayout: invalid P_H component.");
	}
	return 6 * ndofs_ + 3 * ndofs_ + comp * ndofs_;
}

} // namespace maxwell
