[![ubuntu-gnu](https://github.com/OpenSEMBA/dgtd/actions/workflows/ubuntu-gnu.yml/badge.svg)](https://github.com/OpenSEMBA/dgtd/actions/workflows/ubuntu-gnu.yml)

# semba-dgtd

Maxwell curl-equation solver using discontinuous Galerkin methods (OpenSEMBA / UGR).

## Repository layout

| Path | Contents |
|------|----------|
| `src/` | Driver, evolution operators, DG components, MFEM extensions |
| `external/mfem-geg/` | Required MFEM fork (submodule) |
| `testData/maxwellInputs/` | Example simulation cases (JSON + mesh per folder) |
| `test/` | Unit and integration tests (GoogleTest) |
| `docs/` | Input format, tools, and feature design notes |
| `pythonBindings/` | Optional Python bindings |
| `AGENTS.md` | Project guide for AI coding agents |
| `CLAUDE.md` | General behavioral guidelines for LLM-assisted editing |

## Compiling

**Requirements:** CMake ≥ 3.25.2, [vcpkg](https://github.com/microsoft/vcpkg) (`$VCPKG_ROOT`), Ninja, GCC.

vcpkg installs (see [vcpkg.json](vcpkg.json)): `eigen3`, `gtest`, `fftw3`, `nlohmann-json`.

### MPI builds

Requires METIS 5 and HYPRE built from source; set `METIS_DIR` and `HYPRE_DIR` before configuring.

**METIS:**
```sh
wget https://github.com/mfem/tpls/raw/gh-pages/metis-5.1.0.tar.gz
tar -zxvf metis-5.1.0.tar.gz
cd metis-5.1.0
make BUILDDIR=lib config
make BUILDDIR=lib
cp lib/libmetis/libmetis.a lib/
export METIS_DIR=$PWD
```

**HYPRE:**
```sh
wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.31.0.tar.gz
tar -zxvf v2.31.0.tar.gz
cd hypre-2.31.0/src
./configure --disable-fortran
make -j $(nproc)
export HYPRE_DIR=$PWD/hypre
```

### CUDA builds

Requires HYPRE built with CUDA (CMake, not autoconf `./configure`).

**Prerequisites:** [CUDA toolkit](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/); set `-DCMAKE_CUDA_ARCHITECTURES` to your GPU (e.g. `89` for RTX 40-series, `86` for RTX 30-series).

**HYPRE with CUDA:**
```sh
wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.31.0.tar.gz
tar -zxvf v2.31.0.tar.gz
cmake -S hypre-2.31.0/src -B hypre-cuda-build \
      -DHYPRE_WITH_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=<GPU_ARCH> \
      -DCMAKE_INSTALL_PREFIX=$HOME/hypre-cuda-install
cmake --build hypre-cuda-build -j $(nproc)
cmake --install hypre-cuda-build
export HYPRE_CUDA_DIR=$HOME/hypre-cuda-install
```

Then (with `METIS_DIR` set):
```sh
cmake --preset gnu-release-cuda
cmake --build --preset build-gnu-release-cuda --parallel
```

### CMake presets

| Preset | Description |
|--------|-------------|
| `gnu-debug-mpi` | Debug, MPI, OpenMP |
| `gnu-release-mpi` | Release, MPI, OpenMP |
| `gnu-debug-cuda` | Debug, MPI, OpenMP, CUDA (gcc-12) |
| `gnu-release-cuda` | Release, MPI, OpenMP, CUDA (gcc-12) |

```sh
cmake --preset gnu-release-mpi
cmake --build --preset build-gnu-release-mpi --parallel
```

### MFEM

Built automatically from `external/mfem-geg`. For an external install: `-DSEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY=OFF`.

> **Warning:** Use only [OpenSEMBA/mfem-geg](https://github.com/OpenSEMBA/mfem-geg). Upstream MFEM will not build this project.

## Running a case

Cases live under `testData/maxwellInputs/<case_name>/` with matching `<case_name>.json` and mesh file.

Full JSON reference: **[docs/json-input-format.md](docs/json-input-format.md)**.

Example:
```sh
./build/gnu-release-mpi/bin/opensemba_dgtd testData/maxwellInputs/1D_PEC/1D_PEC.json
```

(Confirm binary name/path for your preset.)

Exports appear under `Exports/` by run mode and case name.

## Further documentation

| Topic | Link |
|-------|------|
| All docs | [docs/README.md](docs/README.md) |
| JSON input | [docs/json-input-format.md](docs/json-input-format.md) |
| MOR → ParaView | [docs/mor2paraview.md](docs/mor2paraview.md) |
| Volumetric PML design | [docs/pml/README.md](docs/pml/README.md) |
| AI / agent context | [AGENTS.md](AGENTS.md) |

## Funding

- Spanish Ministry of Science and Innovation (MICIN/AEI) (Grant PID2022-137495OB-C31).
- European Union, HECATE project (HE-HORIZON-JU-Clean-Aviation-2022-01).
- European Union, FEDER 2020 (B-TIC-700-UGR20).
