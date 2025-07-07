#include "GlobalEvolution.h"
#include <mfem.hpp>
#include "general/forall.hpp"

namespace maxwell{

void load_in_to_eh_gpu(mfem::Vector& in, 
                        std::array<mfem::Vector, 3>& eOld, 
                        std::array<mfem::Vector, 3>& hOld, 
                        const int ndofs)
{
    eOld[0].MakeRef(in,  0      * ndofs, ndofs);
    eOld[1].MakeRef(in,  1      * ndofs, ndofs);
    eOld[2].MakeRef(in,  2      * ndofs, ndofs);
    hOld[0].MakeRef(in, (0 + 3) * ndofs, ndofs);
    hOld[1].MakeRef(in, (1 + 3) * ndofs, ndofs);
    hOld[2].MakeRef(in, (2 + 3) * ndofs, ndofs);
}

void load_eh_to_innew_gpu(const std::array<mfem::Vector, 3>& eOld,
                          const std::array<mfem::Vector, 3>& hOld,
                          mfem::Vector& inNew,
                          int ndofs,
                          int nbrSize)
{

    const double* eOld0 = eOld[0].Read();
    const double* eOld1 = eOld[1].Read();
    const double* eOld2 = eOld[2].Read();
    const double* hOld0 = hOld[0].Read();
    const double* hOld1 = hOld[1].Read();
    const double* hOld2 = hOld[2].Read();

    double* inNew_d = inNew.Write();

    const int blockSize = ndofs + nbrSize;

    mfem::forall(ndofs, [=] MFEM_DEVICE(int v) {
        inNew_d[0 * blockSize + v] = eOld0[v];
        inNew_d[1 * blockSize + v] = eOld1[v];
        inNew_d[2 * blockSize + v] = eOld2[v];
        inNew_d[(0 + 3) * blockSize + v] = hOld0[v];
        inNew_d[(1 + 3) * blockSize + v] = hOld1[v];
        inNew_d[(2 + 3) * blockSize + v] = hOld2[v];
    });
    cudaDeviceSynchronize(); // wait for completion
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
        inNew_d[(0 + 3) * blockSize + ndofs + v] = hOldNbr0[v];
        inNew_d[(1 + 3) * blockSize + ndofs + v] = hOldNbr1[v];
        inNew_d[(2 + 3) * blockSize + ndofs + v] = hOldNbr2[v];
    });
}

}