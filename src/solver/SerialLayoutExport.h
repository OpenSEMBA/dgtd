#pragma once

#include <filesystem>
#include <vector>

#include "mfem.hpp"

namespace maxwell {

/// Owned-dof map from a uniform DG space onto the serial mesh numbering.
/// `mfem_to_serial[GetMyDofOffset() + local_dof]` is the serial scalar dof.
struct SerialDGDofMap {
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0;
    int n_ranks = 1;
    int n_local = 0;
    int n_serial = 0;
    int dpe = 0;
    std::vector<int> local_to_serial;
    std::vector<int> mfem_to_serial;
    std::vector<int> rank_counts;
    std::vector<HYPRE_BigInt> rank_offsets;
};

/// `element_partition[serial_element] == rank`, the array stored on `Model`.
SerialDGDofMap buildSerialDGDofMap(
    const mfem::ParFiniteElementSpace& fes,
    const mfem::Array<int>& element_partition);

/// Gather a rank-local `6 N_local` by `6 (N_local + N_nbr)` matrix into the
/// serial square CSR and write it with `PrintCSR2` on rank 0.
void gatherSerialCSR(
    const mfem::SparseMatrix& local,
    const SerialDGDofMap& map,
    const mfem::Array<HYPRE_BigInt>* face_nbr_glob_dof_map,
    int n_face_nbr_dofs,
    const std::filesystem::path& file_path,
    const char* label);

/// Gather the EH block (`6 N_local`, component-major) into serial order.
/// `serial_eh` is filled on rank 0 and cleared on every other rank.
void gatherSerialEHVector(
    const SerialDGDofMap& map,
    const double* local_eh,
    std::vector<double>& serial_eh);

} // namespace maxwell
