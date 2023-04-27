# maxwell-dgtd
Maxwell's curl equations solver using discontinuous Galerkin methods


## Compiling

Compilation needs vcpkg with the following packages:
- eigen3
- gtest
- metis (for MPI)
- hypre (for MPI)

### OpenMP and MPI in Windows

- OpenMP requires an LLVM compiler. It has been tested Intel OneAPI (Base kit and HPC kit). To compile, use the following CMake Command Arguments when compiling MFEM and maxwell dgtd:
    ```
    -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icx
    ```
- MPI requires mfem to be compiled with MFEM_USE_MPI which requires Hypre and METIS.
  - Option 1, install with vcpkg. make sure to mark METIS_VERSION_5 and that th
  - Option 2, compiling from sources: Hypre can be cloned from https://github.com/hypre-space/hypre. It must be compiled **and installed** using a CMakeLists.txt in hypre/src. It is possible that the following variable has to be manually set:
    ```
        HYPRE_LIBRARIES=[PATH TO lib]\HYPRE.lib
    ```
    METIS 5.0.1 can be compiled in Windows with CMake. For VS2022+, comment out the line 
    ```
        #define rint(x) ((idx_t)((x)+0.5))
    ```
    in metis\GKlib\gk_arch.h

### OpenMP and MPI Linux

- Tested to compile with gcc.
- OpenMP has been tested to work.