# maxwell
Maxwell's curl equations solver using discontinuous Galerkin methods


## Compiling

### Windows

- OpenMP requires an LLVM compiler. It has been tested Intel OneAPI HPC.
- MPI requires mfem to be compiled with MFEM_USE_MPI which requires Hypre and METIS.
  - Hypre can be cloned from https://github.com/hypre-space/hypre. It must be compiled **and installed** using a CMakeLists.txt in hypre/src. It is possible that the following variable has to be manually set:
    ```
        HYPRE_LIBRARIES=[PATH TO lib]\HYPRE.lib
    ```
  - METIS 5.0.1 can be compiled in Windows with CMake. For VS2022+, comment out the line 
    ```
        #define rint(x) ((idx_t)((x)+0.5))
    ```
    in metis\GKlib\gk_arch.h

### Linux

- Tested to compile with gcc.
- OpenMP has been tested
- MPI requires