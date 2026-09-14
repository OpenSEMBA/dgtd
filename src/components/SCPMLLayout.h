#pragma once

#include "PMLProperties.h"
#include "Types.h"

#include <stdexcept>
#include <vector>

namespace maxwell {

/// Layout of Bagci/Chen SC-PML auxiliaries in the extended ODE state:
///   [Ex Ey Ez Hx Hy Hz | PEx PEy PEz | PHx PHy PHz]
/// n_aux = 6 * ndofs when any Cartesian PML region exists; else 0.
class SCPMLLayout {
public:
	SCPMLLayout() = default;

	SCPMLLayout(int ndofs, const std::vector<PMLProperties>& regions, int mesh_dim);

	int ndofs() const { return ndofs_; }
	int nAux() const { return n_aux_; }
	bool active() const { return n_aux_ > 0; }

	int pEOffset(Direction comp) const;
	int pHOffset(Direction comp) const;

private:
	int ndofs_ = 0;
	int n_aux_ = 0;
};

inline int computePMLAuxSize(
	const std::vector<PMLProperties>& regions, int ndofs, int mesh_dim)
{
	return SCPMLLayout(ndofs, regions, mesh_dim).nAux();
}

} // namespace maxwell
