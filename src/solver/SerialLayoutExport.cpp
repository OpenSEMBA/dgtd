#include "SerialLayoutExport.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace maxwell {

namespace {

struct Triplet {
    int row;
    int col;
    mfem::real_t val;
};

void require(bool ok, const std::string& message)
{
    if (!ok) {
        throw std::runtime_error(message);
    }
}

int checkIntProduct(int a, int b, const char* what)
{
    require(a >= 0 && b >= 0, std::string(what) + " is negative.");
    if (a == 0 || b == 0) {
        return 0;
    }
    require(a <= std::numeric_limits<int>::max() / b,
            std::string(what) + " does not fit in a 32-bit SparseMatrix index.");
    return a * b;
}

} // namespace

SerialDGDofMap buildSerialDGDofMap(
    const mfem::ParFiniteElementSpace& fes,
    const mfem::Array<int>& element_partition)
{
    SerialDGDofMap map;
    map.comm = fes.GetParMesh()->GetComm();
    map.rank = fes.GetParMesh()->GetMyRank();
    map.n_ranks = fes.GetParMesh()->GetNRanks();
    map.n_local = fes.GetNDofs();

    require(element_partition.Size() > 0,
            "Element partition is empty. MPI MOR export needs Model::elementPartition().");
    require(map.n_ranks > 0, "Parallel mesh has no ranks.");

    const int n_serial_elem = element_partition.Size();
    for (int e = 0; e < n_serial_elem; ++e) {
        const int owner = element_partition[e];
        require(owner >= 0 && owner < map.n_ranks,
                "Element partition rank is outside this job.");
    }

    int local_dpe = 0;
    if (fes.GetNE() > 0) {
        local_dpe = fes.GetFE(0)->GetDof();
        require(local_dpe > 0, "DG element has no degrees of freedom.");
        require(map.n_local == fes.GetNE() * local_dpe,
                "Local DG dof count is not one block per owned element.");
    }
    int dpe_max = 0;
    MPI_Allreduce(&local_dpe, &dpe_max, 1, MPI_INT, MPI_MAX, map.comm);
    const int dpe_for_min = (local_dpe == 0) ? dpe_max : local_dpe;
    int dpe_min = 0;
    MPI_Allreduce(&dpe_for_min, &dpe_min, 1, MPI_INT, MPI_MIN, map.comm);
    require(dpe_max > 0 && dpe_max == dpe_min,
            "DG degrees of freedom per element differ across ranks.");
    map.dpe = dpe_max;

    map.n_serial = checkIntProduct(n_serial_elem, map.dpe, "Serial DG dof count");

    std::vector<int> owned;
    owned.reserve(static_cast<std::size_t>(fes.GetNE()));
    for (int e = 0; e < n_serial_elem; ++e) {
        if (element_partition[e] == map.rank) {
            owned.push_back(e);
        }
    }
    require(static_cast<int>(owned.size()) == fes.GetNE(),
            "Local element count does not match the element partition.");

    map.local_to_serial.assign(static_cast<std::size_t>(map.n_local), -1);
    std::vector<char> seen(static_cast<std::size_t>(map.n_local), 0);
    mfem::Array<int> dofs;
    for (int le = 0; le < fes.GetNE(); ++le) {
        fes.GetElementDofs(le, dofs);
        require(dofs.Size() == map.dpe,
                "Element dof count does not match the DG space.");
        const int serial_e = owned[static_cast<std::size_t>(le)];
        for (int k = 0; k < map.dpe; ++k) {
            const int ldof = mfem::UnsignIndex(dofs[k]);
            require(ldof == le * map.dpe + k,
                    "DG dof numbering is not element-major, so MPI MOR export "
                    "cannot match a one-rank run.");
            require(ldof >= 0 && ldof < map.n_local && !seen[static_cast<std::size_t>(ldof)],
                    "Local DG dof was visited twice.");
            seen[static_cast<std::size_t>(ldof)] = 1;
            map.local_to_serial[static_cast<std::size_t>(ldof)] =
                serial_e * map.dpe + k;
        }
    }

    map.rank_counts.assign(static_cast<std::size_t>(map.n_ranks), 0);
    MPI_Allgather(&map.n_local, 1, MPI_INT,
                  map.rank_counts.data(), 1, MPI_INT, map.comm);

    HYPRE_BigInt my_offset = fes.GetMyDofOffset();
    map.rank_offsets.assign(static_cast<std::size_t>(map.n_ranks), 0);
    MPI_Allgather(&my_offset, 1, HYPRE_MPI_BIG_INT,
                  map.rank_offsets.data(), 1, HYPRE_MPI_BIG_INT, map.comm);

    int gathered = 0;
    std::vector<int> displs(static_cast<std::size_t>(map.n_ranks), 0);
    for (int r = 0; r < map.n_ranks; ++r) {
        require(map.rank_counts[static_cast<std::size_t>(r)] >= 0,
                "Negative local dof count in the serial DG map.");
        displs[static_cast<std::size_t>(r)] = gathered;
        require(gathered <= std::numeric_limits<int>::max()
                    - map.rank_counts[static_cast<std::size_t>(r)],
                "Serial DG dof count does not fit in a 32-bit index.");
        gathered += map.rank_counts[static_cast<std::size_t>(r)];
    }
    require(gathered == map.n_serial,
            "Sum of local DG dofs does not equal the serial dof count.");

    std::vector<int> all_local_to_serial(static_cast<std::size_t>(gathered));
    const int dummy_send = 0;
    const int* send_ptr = map.n_local > 0 ? map.local_to_serial.data() : &dummy_send;
    MPI_Allgatherv(send_ptr, map.n_local, MPI_INT,
                   all_local_to_serial.data(), map.rank_counts.data(), displs.data(),
                   MPI_INT, map.comm);

    map.mfem_to_serial.assign(static_cast<std::size_t>(map.n_serial), -1);
    std::vector<char> slot_filled(static_cast<std::size_t>(map.n_serial), 0);
    std::vector<char> serial_seen(static_cast<std::size_t>(map.n_serial), 0);
    for (int r = 0; r < map.n_ranks; ++r) {
        const int count = map.rank_counts[static_cast<std::size_t>(r)];
        if (count == 0) {
            continue;
        }
        const HYPRE_BigInt off = map.rank_offsets[static_cast<std::size_t>(r)];
        require(off >= 0 && off <= static_cast<HYPRE_BigInt>(map.n_serial - count),
                "MFEM dof offset does not tile the serial DG dof range.");
        for (int i = 0; i < count; ++i) {
            const int slot = static_cast<int>(off) + i;
            require(!slot_filled[static_cast<std::size_t>(slot)],
                    "MFEM dof offsets overlap.");
            slot_filled[static_cast<std::size_t>(slot)] = 1;
            const int serial = all_local_to_serial[
                static_cast<std::size_t>(displs[static_cast<std::size_t>(r)] + i)];
            require(serial >= 0 && serial < map.n_serial &&
                        !serial_seen[static_cast<std::size_t>(serial)],
                    "Serial DG dof map is not a permutation.");
            serial_seen[static_cast<std::size_t>(serial)] = 1;
            map.mfem_to_serial[static_cast<std::size_t>(slot)] = serial;
        }
    }
    for (int i = 0; i < map.n_serial; ++i) {
        require(slot_filled[static_cast<std::size_t>(i)],
                "MFEM dof offsets leave a gap in the serial DG dof range.");
    }
    return map;
}

void gatherSerialCSR(
    const mfem::SparseMatrix& local,
    const SerialDGDofMap& map,
    const mfem::Array<HYPRE_BigInt>* face_nbr_glob_dof_map,
    int n_face_nbr_dofs,
    const std::filesystem::path& file_path,
    const char* label)
{
    static_assert(sizeof(mfem::real_t) == sizeof(double),
                  "MOR export expects double matrix entries.");

    require(local.Finalized(), "Operator export requires a finalized matrix.");
    require(n_face_nbr_dofs >= 0, "Face-neighbor dof count is not initialized.");
    const int n_local = map.n_local;
    const int stride = n_local + n_face_nbr_dofs;
    const int height = local.Height();
    const int width = local.Width();
    if (n_local == 0) {
        require(height == 0 && width == 0,
                "Rank with no dofs holds a non-empty operator block.");
    } else {
        require(height == checkIntProduct(6, n_local, "Local operator height"),
                "Local operator height is not 6 N_local.");
        require(width == checkIntProduct(6, stride, "Local operator width"),
                "Local operator width is not 6 (N_local + N_nbr).");
    }
    if (n_face_nbr_dofs > 0) {
        require(face_nbr_glob_dof_map != nullptr
                    && face_nbr_glob_dof_map->Size() >= n_face_nbr_dofs,
                "Face-neighbor global dof map is missing.");
    }

    const int n_rows = checkIntProduct(6, map.n_serial, "Serial operator dimension");
    const int* I = local.HostReadI();
    const int* J = local.HostReadJ();
    const mfem::real_t* A = local.HostReadData();
    const int nnz_local = (height == 0) ? 0 : I[height];
    const HYPRE_BigInt* nbr = (n_face_nbr_dofs > 0)
        ? face_nbr_glob_dof_map->HostRead() : nullptr;

    std::vector<int> send_rows;
    std::vector<int> send_cols;
    std::vector<double> send_vals;
    send_rows.reserve(static_cast<std::size_t>(nnz_local));
    send_cols.reserve(static_cast<std::size_t>(nnz_local));
    send_vals.reserve(static_cast<std::size_t>(nnz_local));

    for (int row = 0; row < height; ++row) {
        const int comp = row / n_local;
        const int ldof = row - comp * n_local;
        const int serial_row =
            comp * map.n_serial + map.local_to_serial[static_cast<std::size_t>(ldof)];
        for (int k = I[row]; k < I[row + 1]; ++k) {
            const int col = J[k];
            require(col >= 0 && col < width, "Local operator column is out of range.");
            const int col_comp = col / stride;
            const int local_col = col - col_comp * stride;
            int serial_scalar = -1;
            if (local_col < n_local) {
                serial_scalar = map.local_to_serial[static_cast<std::size_t>(local_col)];
            } else {
                const int nb = local_col - n_local;
                require(nb >= 0 && nb < n_face_nbr_dofs,
                        "Halo column is outside the face-neighbor block.");
                const HYPRE_BigInt global_dof = nbr[nb];
                require(global_dof >= 0
                            && global_dof < static_cast<HYPRE_BigInt>(map.n_serial),
                        "Face-neighbor dof is outside the serial DG range.");
                serial_scalar = map.mfem_to_serial[static_cast<std::size_t>(global_dof)];
            }
            require(serial_scalar >= 0 && serial_scalar < map.n_serial,
                    "Serial column dof is out of range.");
            send_rows.push_back(serial_row);
            send_cols.push_back(col_comp * map.n_serial + serial_scalar);
            send_vals.push_back(static_cast<double>(A[k]));
        }
    }

    const int send_nnz = static_cast<int>(send_rows.size());
    std::vector<int> recv_counts(static_cast<std::size_t>(map.n_ranks), 0);
    MPI_Gather(&send_nnz, 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT, 0, map.comm);

    std::vector<int> recv_displs(static_cast<std::size_t>(map.n_ranks), 0);
    int total_nnz = 0;
    if (map.rank == 0) {
        for (int r = 0; r < map.n_ranks; ++r) {
            require(recv_counts[static_cast<std::size_t>(r)] >= 0,
                    "Negative nonzero count while gathering an operator.");
            recv_displs[static_cast<std::size_t>(r)] = total_nnz;
            require(total_nnz <= std::numeric_limits<int>::max()
                        - recv_counts[static_cast<std::size_t>(r)],
                    "Gathered operator has too many nonzeros for a 32-bit CSR.");
            total_nnz += recv_counts[static_cast<std::size_t>(r)];
        }
    }

    std::vector<int> all_rows(static_cast<std::size_t>(map.rank == 0 ? std::max(total_nnz, 1) : 1));
    std::vector<int> all_cols(static_cast<std::size_t>(map.rank == 0 ? std::max(total_nnz, 1) : 1));
    std::vector<double> all_vals(static_cast<std::size_t>(map.rank == 0 ? std::max(total_nnz, 1) : 1));
    const int dummy_int = 0;
    const double dummy_real = 0.0;
    MPI_Gatherv(send_nnz > 0 ? send_rows.data() : &dummy_int, send_nnz, MPI_INT,
                all_rows.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                0, map.comm);
    MPI_Gatherv(send_nnz > 0 ? send_cols.data() : &dummy_int, send_nnz, MPI_INT,
                all_cols.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                0, map.comm);
    MPI_Gatherv(send_nnz > 0 ? send_vals.data() : &dummy_real, send_nnz, MPI_DOUBLE,
                all_vals.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                0, map.comm);

    int write_ok = 1;
    if (map.rank == 0) {
        try {
            std::vector<Triplet> triplets(static_cast<std::size_t>(total_nnz));
            for (int i = 0; i < total_nnz; ++i) {
                triplets[static_cast<std::size_t>(i)] = Triplet{
                    all_rows[static_cast<std::size_t>(i)],
                    all_cols[static_cast<std::size_t>(i)],
                    all_vals[static_cast<std::size_t>(i)]};
                require(triplets[static_cast<std::size_t>(i)].row >= 0
                            && triplets[static_cast<std::size_t>(i)].row < n_rows
                            && triplets[static_cast<std::size_t>(i)].col >= 0
                            && triplets[static_cast<std::size_t>(i)].col < n_rows,
                        "Gathered operator entry is outside the serial matrix.");
            }
            std::sort(triplets.begin(), triplets.end(),
                      [](const Triplet& a, const Triplet& b) {
                          if (a.row != b.row) {
                              return a.row < b.row;
                          }
                          return a.col < b.col;
                      });

            std::vector<int> row_nnz(static_cast<std::size_t>(n_rows), 0);
            std::vector<int> cols;
            std::vector<mfem::real_t> vals;
            cols.reserve(triplets.size());
            vals.reserve(triplets.size());
            for (std::size_t i = 0; i < triplets.size();) {
                const int row = triplets[i].row;
                const int col = triplets[i].col;
                mfem::real_t val = triplets[i].val;
                std::size_t j = i + 1;
                while (j < triplets.size()
                       && triplets[j].row == row
                       && triplets[j].col == col) {
                    val += triplets[j].val;
                    ++j;
                }
                if (val != 0.0) {
                    cols.push_back(col);
                    vals.push_back(val);
                    ++row_nnz[static_cast<std::size_t>(row)];
                }
                i = j;
            }

            const int kept = static_cast<int>(cols.size());
            std::unique_ptr<mfem::SparseMatrix> serial;
            if (kept == 0) {
                serial = std::make_unique<mfem::SparseMatrix>(n_rows, n_rows);
                serial->Finalize();
            } else {
                int* C_i = mfem::Memory<int>(n_rows + 1);
                int* C_j = mfem::Memory<int>(kept);
                mfem::real_t* C_a = mfem::Memory<mfem::real_t>(kept);
                C_i[0] = 0;
                for (int row = 0; row < n_rows; ++row) {
                    C_i[row + 1] = C_i[row] + row_nnz[static_cast<std::size_t>(row)];
                }
                std::copy(cols.begin(), cols.end(), C_j);
                std::copy(vals.begin(), vals.end(), C_a);
                serial.reset(new mfem::SparseMatrix(C_i, C_j, C_a, n_rows, n_rows));
            }

            std::filesystem::create_directories(file_path.parent_path());
            std::ofstream ofs(file_path);
            if (!ofs.is_open()) {
                throw std::runtime_error(
                    "Could not open file for writing: " + file_path.string());
            }
            serial->PrintCSR2(ofs);
            std::cout << label << " exported to " << file_path << std::endl;
        } catch (const std::exception& ex) {
            write_ok = 0;
            std::cerr << ex.what() << std::endl;
        }
    }

    MPI_Bcast(&write_ok, 1, MPI_INT, 0, map.comm);
    require(write_ok != 0,
            "Failed to write gathered operator: " + file_path.string());
}

void gatherSerialEHVector(
    const SerialDGDofMap& map,
    const double* local_eh,
    std::vector<double>& serial_eh)
{
    const int send_count = checkIntProduct(6, map.n_local, "Local EH vector");
    const int serial_count = checkIntProduct(6, map.n_serial, "Serial EH vector");

    std::vector<int> recv_counts(static_cast<std::size_t>(map.n_ranks), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(map.n_ranks), 0);
    int total = 0;
    for (int r = 0; r < map.n_ranks; ++r) {
        const int count = checkIntProduct(6, map.rank_counts[static_cast<std::size_t>(r)],
                                          "Per-rank EH vector");
        recv_counts[static_cast<std::size_t>(r)] = count;
        recv_displs[static_cast<std::size_t>(r)] = total;
        require(total <= std::numeric_limits<int>::max() - count,
                "Serial EH vector does not fit in a 32-bit length.");
        total += count;
    }
    require(total == serial_count, "Gathered EH length is not 6 N_serial.");

    std::vector<double> recv(static_cast<std::size_t>(
        map.rank == 0 ? std::max(total, 1) : 1));
    const double dummy = 0.0;
    const double* send_ptr = send_count > 0 ? local_eh : &dummy;
    require(send_count == 0 || local_eh != nullptr, "Local EH pointer is null.");
    MPI_Gatherv(send_ptr, send_count, MPI_DOUBLE,
                recv.data(), recv_counts.data(), recv_displs.data(), MPI_DOUBLE,
                0, map.comm);

    if (map.rank != 0) {
        serial_eh.clear();
        return;
    }

    serial_eh.assign(static_cast<std::size_t>(serial_count), 0.0);
    for (int r = 0; r < map.n_ranks; ++r) {
        const int n = map.rank_counts[static_cast<std::size_t>(r)];
        const int off = static_cast<int>(map.rank_offsets[static_cast<std::size_t>(r)]);
        const double* block = recv.data() + recv_displs[static_cast<std::size_t>(r)];
        for (int comp = 0; comp < 6; ++comp) {
            for (int ldof = 0; ldof < n; ++ldof) {
                const int serial = map.mfem_to_serial[static_cast<std::size_t>(off + ldof)];
                serial_eh[static_cast<std::size_t>(comp * map.n_serial + serial)] =
                    block[comp * n + ldof];
            }
        }
    }
}

} // namespace maxwell
