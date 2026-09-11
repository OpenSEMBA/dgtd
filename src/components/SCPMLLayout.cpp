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

	bool any_pml = false;
	for (const auto& props : regions) {
		if (props.uniaxial_radial) {
			if (mesh_dim < 2) {
				throw std::runtime_error(
					"SCPMLLayout: uniaxial radial PML requires dim >= 2.");
			}
			any_pml = true;
			continue;
		}
		for (Direction d : props.active_axes) {
			if (d < 0 || d >= mesh_dim) {
				throw std::runtime_error(
					"SCPMLLayout: active_axes exceeds mesh dimension.");
			}
			any_pml = true;
		}
	}
	if (any_pml) {
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
