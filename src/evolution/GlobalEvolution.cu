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

}