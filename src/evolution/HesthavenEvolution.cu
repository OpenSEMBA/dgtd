#include "HesthavenEvolution.h"
#include <mfem.hpp>
#include "general/forall.hpp"

namespace maxwell {

void hesthaven_compute_jumps_gpu(HesthavenGPUData& gpu,
                                 const mfem::Vector& in,
                                 const std::array<mfem::ParGridFunction, 3>& eOld,
                                 const std::array<mfem::ParGridFunction, 3>& hOld,
                                 int ndofs)
{
    const int jumps_size = gpu.jumps_size;
    const double* in_d = in.Read();
    double* jumps_e = gpu.d_jumps_e.Write();
    double* jumps_h = gpu.d_jumps_h.Write();

    const int* jump_minus = gpu.d_jump_minus.Read();
    const int* jump_plus = gpu.d_jump_plus.Read();
    const uint8_t* plus_is_nbr = gpu.d_jump_plus_is_nbr.Read();
    const bool use_nbr_flags = (gpu.d_jump_plus_is_nbr.Size() == jumps_size);

    const double* e_nbr0 = eOld[0].FaceNbrData().Read();
    const double* e_nbr1 = eOld[1].FaceNbrData().Read();
    const double* e_nbr2 = eOld[2].FaceNbrData().Read();
    const double* h_nbr0 = hOld[0].FaceNbrData().Read();
    const double* h_nbr1 = hOld[1].FaceNbrData().Read();
    const double* h_nbr2 = hOld[2].FaceNbrData().Read();

    mfem::forall(jumps_size, [=] MFEM_DEVICE(int v) {
        const int minus = jump_minus[v];
        const int plus = jump_plus[v];
        const bool from_nbr = use_nbr_flags && (plus_is_nbr[v] != 0);

        const double e_minus = in_d[0 * ndofs + minus];
        const double e_plus = from_nbr ? e_nbr0[plus] : in_d[0 * ndofs + plus];
        jumps_e[0 * jumps_size + v] = e_plus - e_minus;

        const double ey_minus = in_d[1 * ndofs + minus];
        const double ey_plus = from_nbr ? e_nbr1[plus] : in_d[1 * ndofs + plus];
        jumps_e[1 * jumps_size + v] = ey_plus - ey_minus;

        const double ez_minus = in_d[2 * ndofs + minus];
        const double ez_plus = from_nbr ? e_nbr2[plus] : in_d[2 * ndofs + plus];
        jumps_e[2 * jumps_size + v] = ez_plus - ez_minus;

        const double h_minus = in_d[3 * ndofs + minus];
        const double h_plus = from_nbr ? h_nbr0[plus] : in_d[3 * ndofs + plus];
        jumps_h[0 * jumps_size + v] = h_plus - h_minus;

        const double hy_minus = in_d[4 * ndofs + minus];
        const double hy_plus = from_nbr ? h_nbr1[plus] : in_d[4 * ndofs + plus];
        jumps_h[1 * jumps_size + v] = hy_plus - hy_minus;

        const double hz_minus = in_d[5 * ndofs + minus];
        const double hz_plus = from_nbr ? h_nbr2[plus] : in_d[5 * ndofs + plus];
        jumps_h[2 * jumps_size + v] = hz_plus - hz_minus;
    });
}

void hesthaven_apply_bc_gpu(HesthavenGPUData& gpu,
                            const mfem::Vector& in,
                            int ndofs)
{
    const int jumps_size = gpu.jumps_size;
    const double* in_d = in.Read();
    double* jumps_e = gpu.d_jumps_e.Write();
    double* jumps_h = gpu.d_jumps_h.Write();

    const int n_true = gpu.n_bc_true;
    const int* jump_out = gpu.d_bc_true_jump_out.Read();
    const int* vmap_in = gpu.d_bc_true_vmap_in.Read();
    const int* comp = gpu.d_bc_true_comp.Read();
    const double* e_coeff = gpu.d_bc_true_e_coeff.Read();
    const double* h_coeff = gpu.d_bc_true_h_coeff.Read();

    mfem::forall(n_true, [=] MFEM_DEVICE(int i) {
        const int c = comp[i];
        const int vin = vmap_in[i];
        const int jout = jump_out[i];
        jumps_e[c * jumps_size + jout] = e_coeff[i] * in_d[c * ndofs + vin];
        jumps_h[c * jumps_size + jout] = h_coeff[i] * in_d[(3 + c) * ndofs + vin];
    });

    const int n_int = gpu.n_bc_int;
    const int* jump_out1 = gpu.d_bc_int_jump_out1.Read();
    const int* jump_out2 = gpu.d_bc_int_jump_out2.Read();
    const int* vmap_in1 = gpu.d_bc_int_vmap_in1.Read();
    const int* vmap_in2 = gpu.d_bc_int_vmap_in2.Read();
    const int* comp_i = gpu.d_bc_int_comp.Read();
    const double* e_coeff_i = gpu.d_bc_int_e_coeff.Read();
    const double* h_coeff_i = gpu.d_bc_int_h_coeff.Read();

    mfem::forall(n_int, [=] MFEM_DEVICE(int i) {
        const int c = comp_i[i];
        const double e_val1 = e_coeff_i[i] * in_d[c * ndofs + vmap_in1[i]];
        const double h_val1 = h_coeff_i[i] * in_d[(3 + c) * ndofs + vmap_in1[i]];
        const double e_val2 = e_coeff_i[i] * in_d[c * ndofs + vmap_in2[i]];
        const double h_val2 = h_coeff_i[i] * in_d[(3 + c) * ndofs + vmap_in2[i]];
        jumps_e[c * jumps_size + jump_out1[i]] = e_val1;
        jumps_h[c * jumps_size + jump_out1[i]] = h_val1;
        jumps_e[c * jumps_size + jump_out2[i]] = e_val2;
        jumps_h[c * jumps_size + jump_out2[i]] = h_val2;
    });
}

void hesthaven_scatter_tfsf_to_jumps_gpu(HesthavenGPUData& gpu,
                                         const FieldGridFuncs& tfsf_fields,
                                         int jumps_size)
{
    const int n = gpu.n_tfsf;
    if (n == 0) {
        return;
    }

    const int* jump_sf = gpu.d_tfsf_jump_sf.Read();
    const int* jump_tf = gpu.d_tfsf_jump_tf.Read();
    const int* src_dof = gpu.d_tfsf_src_dof.Read();

    const double* e0 = tfsf_fields[E][X].Read();
    const double* e1 = tfsf_fields[E][Y].Read();
    const double* e2 = tfsf_fields[E][Z].Read();
    const double* h0 = tfsf_fields[H][X].Read();
    const double* h1 = tfsf_fields[H][Y].Read();
    const double* h2 = tfsf_fields[H][Z].Read();

    double* jumps_e = gpu.d_jumps_e.Write();
    double* jumps_h = gpu.d_jumps_h.Write();

    mfem::forall(n, [=] MFEM_DEVICE(int i) {
        const int src = src_dof[i];
        const int jsf = jump_sf[i];
        const int jtf = jump_tf[i];

        const double ee0 = e0[src];
        const double ee1 = e1[src];
        const double ee2 = e2[src];
        const double hh0 = h0[src];
        const double hh1 = h1[src];
        const double hh2 = h2[src];

        jumps_e[0 * jumps_size + jsf] += ee0;
        jumps_e[1 * jumps_size + jsf] += ee1;
        jumps_e[2 * jumps_size + jsf] += ee2;
        jumps_h[0 * jumps_size + jsf] += hh0;
        jumps_h[1 * jumps_size + jsf] += hh1;
        jumps_h[2 * jumps_size + jsf] += hh2;

        jumps_e[0 * jumps_size + jtf] -= ee0;
        jumps_e[1 * jumps_size + jtf] -= ee1;
        jumps_e[2 * jumps_size + jtf] -= ee2;
        jumps_h[0 * jumps_size + jtf] -= hh0;
        jumps_h[1 * jumps_size + jtf] -= hh1;
        jumps_h[2 * jumps_size + jtf] -= hh2;
    });
}

namespace {

// Eigen::MatrixXd is column-major; element (i, j) is at mat[i + j * rows].
MFEM_HOST_DEVICE inline void dense_matvec(int rows, int cols,
                                          const double* mat,
                                          const double* x,
                                          double* y)
{
    for (int i = 0; i < rows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < cols; ++j) {
            sum += mat[i + j * rows] * x[j];
        }
        y[i] = sum;
    }
}

} // namespace

void hesthaven_mult_gpu(HesthavenGPUData& gpu,
                        double alpha,
                        const mfem::Vector& in,
                        const std::array<mfem::ParGridFunction, 3>&,
                        const std::array<mfem::ParGridFunction, 3>&,
                        mfem::Vector& out,
                        int ndofs)
{
    const int n_elem = gpu.n_elements;
    const int jumps_size = gpu.jumps_size;
    const int ndof_el = gpu.ndof_el;
    const int flux_cap = gpu.max_flux_size;
    const int lift_rows = gpu.lift_rows;
    const int lift_cols = gpu.lift_cols;
    const int ws_stride = gpu.workspace_stride;

    const double* in_d = in.Read();
    double* out_d = out.Write();

    const double* jumps_e = gpu.d_jumps_e.Read();
    const double* jumps_h = gpu.d_jumps_h.Read();

    const int* elem_ids = gpu.d_elem_ids.Read();
    const int* elem_dofs = gpu.d_elem_dofs.Read();
    const int* dir_matrix_id = gpu.d_dir_matrix_id.Read();
    const int* flux_size_arr = gpu.d_flux_size.Read();
    const double* normals = gpu.d_normals.Read();
    const double* fscale = gpu.d_fscale.Read();
    const double* matrices = gpu.d_matrices.Read();
    const int* matrix_offsets = gpu.d_matrix_offsets.Read();
    const int* matrix_rows = gpu.d_matrix_rows.Read();
    const int* matrix_cols = gpu.d_matrix_cols.Read();
    const double* ref_lift = gpu.d_ref_lift.Read();
    double* workspace = gpu.d_workspace.Write();

    mfem::forall(n_elem, [=] MFEM_DEVICE(int e) {
        const int elem_id = elem_ids[e];
        const int flux_size = flux_size_arr[e];
        const int jump_base = elem_id * flux_size;
        const int dof_base = e * ndof_el;

        double* ws = workspace + e * ws_stride;
        double* ndotdH = ws;
        double* ndotdE = ws + flux_cap;
        double* je_x = ws + 2 * flux_cap;
        double* je_y = ws + 3 * flux_cap;
        double* je_z = ws + 4 * flux_cap;
        double* jh_x = ws + 5 * flux_cap;
        double* jh_y = ws + 6 * flux_cap;
        double* jh_z = ws + 7 * flux_cap;
        const int field_base = 8 * flux_cap;
        double* field_ex = ws + field_base;
        double* field_ey = ws + field_base + ndof_el;
        double* field_ez = ws + field_base + 2 * ndof_el;
        double* field_hx = ws + field_base + 3 * ndof_el;
        double* field_hy = ws + field_base + 4 * ndof_el;
        double* field_hz = ws + field_base + 5 * ndof_el;
        const int flux_base = field_base + 6 * ndof_el;
        double* flux_hx = ws + flux_base;
        double* flux_hy = ws + flux_base + flux_cap;
        double* flux_hz = ws + flux_base + 2 * flux_cap;
        double* flux_ex = ws + flux_base + 3 * flux_cap;
        double* flux_ey = ws + flux_base + 4 * flux_cap;
        double* flux_ez = ws + flux_base + 5 * flux_cap;
        const int lift_base = flux_base + 6 * flux_cap;
        double* lift_in = ws + lift_base;
        double* lift_out = ws + lift_base + flux_cap;
        double* tmp_y = ws + lift_base + flux_cap + ndof_el;
        double* tmp_z = ws + lift_base + 2 * flux_cap + ndof_el;

        for (int i = 0; i < flux_size; ++i) {
            const int jidx = jump_base + i;
            je_x[i] = jumps_e[0 * jumps_size + jidx];
            je_y[i] = jumps_e[1 * jumps_size + jidx];
            je_z[i] = jumps_e[2 * jumps_size + jidx];
            jh_x[i] = jumps_h[0 * jumps_size + jidx];
            jh_y[i] = jumps_h[1 * jumps_size + jidx];
            jh_z[i] = jumps_h[2 * jumps_size + jidx];
        }

        for (int i = 0; i < flux_size; ++i) {
            ndotdH[i] = 0.0;
            ndotdE[i] = 0.0;
            for (int d = 0; d < 3; ++d) {
                const double nor = normals[e * 3 * flux_cap + d * flux_cap + i];
                ndotdH[i] += nor * (d == 0 ? jh_x[i] : (d == 1 ? jh_y[i] : jh_z[i]));
                ndotdE[i] += nor * (d == 0 ? je_x[i] : (d == 1 ? je_y[i] : je_z[i]));
            }
        }

        for (int i = 0; i < flux_size; ++i) {
            const double norx = normals[e * 3 * flux_cap + 0 * flux_cap + i];
            const double nory = normals[e * 3 * flux_cap + 1 * flux_cap + i];
            const double norz = normals[e * 3 * flux_cap + 2 * flux_cap + i];
            flux_hx[i] = -nory * je_z[i] + norz * je_y[i]
                + alpha * (jh_x[i] - ndotdH[i] * norx);
            flux_hy[i] = -norz * je_x[i] + norx * je_z[i]
                + alpha * (jh_y[i] - ndotdH[i] * nory);
            flux_hz[i] = -norx * je_y[i] + nory * je_x[i]
                + alpha * (jh_z[i] - ndotdH[i] * norz);
            flux_ex[i] = nory * jh_z[i] - norz * jh_y[i]
                + alpha * (je_x[i] - ndotdE[i] * norx);
            flux_ey[i] = norz * jh_x[i] - norx * jh_z[i]
                + alpha * (je_y[i] - ndotdE[i] * nory);
            flux_ez[i] = norx * jh_y[i] - nory * jh_x[i]
                + alpha * (je_z[i] - ndotdE[i] * norz);
        }

        for (int i = 0; i < ndof_el; ++i) {
            const int gi = elem_dofs[dof_base + i];
            field_ex[i] = in_d[0 * ndofs + gi];
            field_ey[i] = in_d[1 * ndofs + gi];
            field_ez[i] = in_d[2 * ndofs + gi];
            field_hx[i] = in_d[3 * ndofs + gi];
            field_hy[i] = in_d[4 * ndofs + gi];
            field_hz[i] = in_d[5 * ndofs + gi];
        }

        for (int x = 0; x < 3; ++x) {
            const int y = (x + 1) % 3;
            const int z = (x + 2) % 3;

            const int mat_y = dir_matrix_id[e * 3 + y];
            const int mat_z = dir_matrix_id[e * 3 + z];
            const double* dir_y = matrices + matrix_offsets[mat_y];
            const double* dir_z = matrices + matrix_offsets[mat_z];

            const double* field_e_ycomp = (y == 0) ? field_ex : (y == 1 ? field_ey : field_ez);
            const double* field_e_zcomp = (z == 0) ? field_ex : (z == 1 ? field_ey : field_ez);
            const double* field_h_ycomp = (y == 0) ? field_hx : (y == 1 ? field_hy : field_hz);
            const double* field_h_zcomp = (z == 0) ? field_hx : (z == 1 ? field_hy : field_hz);
            double* flux_h_xcomp = (x == 0) ? flux_hx : (x == 1 ? flux_hy : flux_hz);
            double* flux_e_xcomp = (x == 0) ? flux_ex : (x == 1 ? flux_ey : flux_ez);

            dense_matvec(matrix_rows[mat_y], matrix_cols[mat_y], dir_y, field_e_zcomp, tmp_y);
            dense_matvec(matrix_rows[mat_z], matrix_cols[mat_z], dir_z, field_e_ycomp, tmp_z);
            for (int i = 0; i < ndof_el; ++i) {
                const int gi = elem_dofs[dof_base + i];
                out_d[(3 + x) * ndofs + gi] += -tmp_y[i] + tmp_z[i];
            }

            for (int i = 0; i < flux_size; ++i) {
                lift_in[i] = flux_h_xcomp[i] * fscale[e * flux_cap + i] * 0.5;
            }
            dense_matvec(lift_rows, lift_cols, ref_lift, lift_in, lift_out);
            for (int i = 0; i < ndof_el; ++i) {
                const int gi = elem_dofs[dof_base + i];
                out_d[(3 + x) * ndofs + gi] += lift_out[i];
            }

            dense_matvec(matrix_rows[mat_y], matrix_cols[mat_y], dir_y, field_h_zcomp, tmp_y);
            dense_matvec(matrix_rows[mat_z], matrix_cols[mat_z], dir_z, field_h_ycomp, tmp_z);
            for (int i = 0; i < ndof_el; ++i) {
                const int gi = elem_dofs[dof_base + i];
                out_d[x * ndofs + gi] += tmp_y[i] - tmp_z[i];
            }

            for (int i = 0; i < flux_size; ++i) {
                lift_in[i] = flux_e_xcomp[i] * fscale[e * flux_cap + i] * 0.5;
            }
            dense_matvec(lift_rows, lift_cols, ref_lift, lift_in, lift_out);
            for (int i = 0; i < ndof_el; ++i) {
                const int gi = elem_dofs[dof_base + i];
                out_d[x * ndofs + gi] += lift_out[i];
            }
        }
    });

    MFEM_STREAM_SYNC;
}

} // namespace maxwell
