#pragma once

#include "PMLProperties.h"
#include "Types.h"

#include <vector>

namespace maxwell {

class Model;

/// Describes stretch-direction ordering and offsets for PML auxiliary DOFs.
/// Layout per stretch slot: [ψ^E_{d,X}, ψ^E_{d,Y}, ψ^E_{d,Z}, ψ^H_{d,X}, ψ^H_{d,Y}, ψ^H_{d,Z}]
/// (each block has ndofs scalar L2 DOFs).
class PMLAuxLayout {
public:
	PMLAuxLayout() = default;
	PMLAuxLayout(const mfem::ParFiniteElementSpace& fes,
	             const std::vector<PMLProperties>& pml_props);

	int ndofs() const { return ndofs_; }
	int fieldBlockSize() const { return field_block_size_; }
	int nAux() const { return n_aux_; }
	int totalSize() const { return field_block_size_ + n_aux_; }
	int auxBlockSize() const { return 6 * ndofs_; }

	int numStretchDirections() const { return static_cast<int>(stretch_dirs_.size()); }
	Direction stretchDirection(int slot) const { return stretch_dirs_[slot]; }

	/// Row offset for ψ^E_{stretch_d, h_comp} (memory in Ĥ_{h_comp}).
	int psiEOffset(Direction stretch_d, Direction h_comp) const;
	/// Row offset for ψ^H_{stretch_d, e_comp} (memory in Ė_{e_comp}).
	int psiHOffset(Direction stretch_d, Direction e_comp) const;

	bool isActiveDirection(Direction d) const;
	int stretchSlot(Direction d) const;

private:
	int ndofs_ = 0;
	int field_block_size_ = 0;
	int n_aux_ = 0;
	std::vector<Direction> stretch_dirs_;
	std::vector<int> dir_to_slot_;
};

int computePMLAuxSize(const Model& model, int ndofs, int mesh_dim);

} // namespace maxwell
