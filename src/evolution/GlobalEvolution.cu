#include "GlobalEvolution.h"
#include <mfem.hpp>
#include "general/forall.hpp"

namespace maxwell{

void load_in_to_eh_gpu(const mfem::Vector& in, 
                        std::array<ParGridFunction, 3>& e,
                        std::array<ParGridFunction, 3>& h,
                        int ndofs)
{
    const double* in_d = in.Read();
    double* ex_d = e[0].Write();
    double* ey_d = e[1].Write();
    double* ez_d = e[2].Write();
    double* hx_d = h[0].Write();
    double* hy_d = h[1].Write();
    double* hz_d = h[2].Write();

    mfem::forall(ndofs, [=] MFEM_DEVICE(int v) {
        ex_d[v] = in_d[0 * ndofs + v];
        ey_d[v] = in_d[1 * ndofs + v];
        ez_d[v] = in_d[2 * ndofs + v];
        hx_d[v] = in_d[3 * ndofs + v];
        hy_d[v] = in_d[4 * ndofs + v];
        hz_d[v] = in_d[5 * ndofs + v];
    });
    cudaDeviceSynchronize();
}

void load_eh_to_innew_gpu(const mfem::Vector& in,
                          mfem::Vector& inNew,
                          int ndofs,
                          int nbrSize)
{

    const double* in_d = in.Read();
    double* inNew_d = inNew.Write();

    const int blockSize = ndofs + nbrSize;

    mfem::forall(ndofs, [=] MFEM_DEVICE(int v) {
        inNew_d[0 * blockSize + v] = in_d[0 * ndofs + v];
        inNew_d[1 * blockSize + v] = in_d[1 * ndofs + v];
        inNew_d[2 * blockSize + v] = in_d[2 * ndofs + v];
        inNew_d[3 * blockSize + v] = in_d[3 * ndofs + v];
        inNew_d[4 * blockSize + v] = in_d[4 * ndofs + v];
        inNew_d[5 * blockSize + v] = in_d[5 * ndofs + v];
    });
    cudaDeviceSynchronize();
}

void load_nbr_to_innew_gpu(const std::array<mfem::ParGridFunction, 3>& eOldNbr,
                   const std::array<mfem::ParGridFunction, 3>& hOldNbr,
                   mfem::Vector& inNew,
                   const int ndofs,
                   const int nbrSize)
{
    const double* eOldNbr0 = eOldNbr[0].FaceNbrData().Read();
    const double* eOldNbr1 = eOldNbr[1].FaceNbrData().Read();
    const double* eOldNbr2 = eOldNbr[2].FaceNbrData().Read();
    const double* hOldNbr0 = hOldNbr[0].FaceNbrData().Read();
    const double* hOldNbr1 = hOldNbr[1].FaceNbrData().Read();
    const double* hOldNbr2 = hOldNbr[2].FaceNbrData().Read();

    double* inNew_d = inNew.Write();

    const int blockSize = ndofs + nbrSize;

    mfem::forall(nbrSize, [=] MFEM_DEVICE(int v) {
        inNew_d[0 * blockSize + ndofs + v] = eOldNbr0[v];
        inNew_d[1 * blockSize + ndofs + v] = eOldNbr1[v];
        inNew_d[2 * blockSize + ndofs + v] = eOldNbr2[v];
        inNew_d[3 * blockSize + ndofs + v] = hOldNbr0[v];
        inNew_d[4 * blockSize + ndofs + v] = hOldNbr1[v];
        inNew_d[5 * blockSize + ndofs + v] = hOldNbr2[v];
    });
    cudaDeviceSynchronize();
}

mfem::Vector load_tfsf_into_single_vector_gpu(const FieldGridFuncs& func)
{
    const int n_fields = 2;
    const int n_dirs = 3;
    const int v_size = func[0][0].Size();
    const int total = n_fields * n_dirs * v_size;

    mfem::Vector res(total);
    res.UseDevice(true);

    auto ex_x = func[0][0].Read();
    auto ex_y = func[0][1].Read();
    auto ex_z = func[0][2].Read();
    auto hx_x = func[1][0].Read();
    auto hx_y = func[1][1].Read();
    auto hx_z = func[1][2].Read();

    auto res_d = res.Write();

    mfem::forall(total, [=] MFEM_HOST_DEVICE (int i) {
        const int v = i % v_size;
        const int d = (i / v_size) % n_dirs;
        const int f = i / (n_dirs * v_size);

        const double val =
            (f == 0 && d == 0) ? ex_x[v] :
            (f == 0 && d == 1) ? ex_y[v] :
            (f == 0 && d == 2) ? ex_z[v] :
            (f == 1 && d == 0) ? hx_x[v] :
            (f == 1 && d == 1) ? hx_y[v] :
                                 hx_z[v];

        res_d[i] = val;
    });

    return res;
}

void load_tfsf_into_out_vector_gpu(const mfem::Array<int>& tfsf_sub_to_parent_ids_, 
                                   const mfem::Vector& tempTFSF,                               
                                         mfem::Vector& out,
                                   const int global_ndofs,
                                   const int tfsf_ndofs)
{
    const int n_fields = 2;
    const int n_dirs = 3;
    const int v_size = tfsf_sub_to_parent_ids_.Size();
    const int total = n_fields * n_dirs * v_size;

    auto tempTFSF_d = tempTFSF.Read();
    auto out_d = out.ReadWrite();
    auto sub_to_parent = tfsf_sub_to_parent_ids_.Read();

    mfem::forall(total, [=] MFEM_HOST_DEVICE (int i) {
        const int v = i % v_size;
        const int d = (i / v_size) % n_dirs;
        const int f = i / (n_dirs * v_size);

        const int outIdx = (f * n_dirs + d) * global_ndofs + sub_to_parent[v];
        const int tempIdx = (f * n_dirs + d) * tfsf_ndofs + v;

        if (tempTFSF_d[tempIdx] != 0.0) {
            out_d[outIdx] -= tempTFSF_d[tempIdx];
        }
    });
}

FieldGridFuncs eval_time_var_field_gpu(const Time time, SourcesManager& sm)
{

    FieldGridFuncs res = sm.evalTimeVarField(time, sm.getGlobalTFSFSpace());
    FieldGridFuncs func_g_sf = res;

    sm.markDoFSforTForSF(res, true);

    if (sm.getTFSFSubMesher().getSFSubMesh() != nullptr)
    {
        sm.markDoFSforTForSF(func_g_sf, false);
        for (int f : {E, H})
        {
            for (int d = 0; d <= Z; ++d)
            {
                mfem::Vector &res_v = res[f][d];
                const mfem::Vector &func_g_v = func_g_sf[f][d];
                const int size = res_v.Size();

                auto res_d = res_v.ReadWrite();
                auto func_g_d = func_g_v.Read();

                mfem::forall(size, [=] MFEM_HOST_DEVICE (int i) {
                    res_d[i] = 0.5 * (res_d[i] - func_g_d[i]);
                });
            }
        }
    }

    return res;
}

}