#pragma once

#include "PMLProperties.h"
#include "Types.h"

#include <array>
#include <stdexcept>
#include <vector>

namespace maxwell {

/// Layout of classical ADE-PML auxiliaries in the extended ODE state:
///   [Ex Ey Ez Hx Hy Hz | J_d0 M_d0 J_d1 M_d1 ... ]
/// with one (J,M) pair per stretch direction in the union of region active_axes.
/// v1 supports uniaxial regions only (one axis per material block).
class ClassicalPMLLayout {
public:
	ClassicalPMLLayout() = default;

	ClassicalPMLLayout(int ndofs, const std::vector<PMLProperties>& regions, int mesh_dim);

	int ndofs() const { return ndofs_; }
	int nAux() const { return n_aux_; }
	int numStretchDirections() const { return static_cast<int>(stretch_dirs_.size()); }
	const std::vector<Direction>& stretchDirections() const { return stretch_dirs_; }

	bool isActiveDirection(Direction d) const;
	int stretchSlot(Direction d) const;

	/// Absolute offsets into the extended state vector (length 6*ndofs + n_aux).
	int jOffset(Direction stretch_d) const;
	int mOffset(Direction stretch_d) const;

private:
	int ndofs_ = 0;
	int n_aux_ = 0;
	std::vector<Direction> stretch_dirs_;
	std::array<int, 3> dir_to_slot_{{-1, -1, -1}};
};

inline int computeClassicalPMLAuxSize(
	const std::vector<PMLProperties>& regions, int ndofs, int mesh_dim)
{
	return ClassicalPMLLayout(ndofs, regions, mesh_dim).nAux();
}

} // namespace maxwell
