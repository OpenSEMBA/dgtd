[![ubuntu-gnu](https://github.com/OpenSEMBA/dgtd/actions/workflows/ubuntu-gnu.yml/badge.svg)](https://github.com/OpenSEMBA/dgtd/actions/workflows/ubuntu-gnu.yml)

# semba-dgtd

Maxwell curl-equation solver using discontinuous Galerkin methods (OpenSEMBA / UGR).

## Repository layout

| Path | Contents |
|------|----------|
| `src/` | Driver, evolution operators, DG components, MFEM extensions |
| `external/mfem-geg/` | Required MFEM fork (submodule) |
| `testData/maxwellInputs/` | Example simulation cases (JSON + mesh per folder) |
| `testData/rcsInputs/` | Offline `opensemba_rcs` JSON (frequency/angle sweeps) |
| `test/` | Unit and integration tests (GoogleTest) |
| `docs/` | Input format, tools, and feature usage notes |
| `pythonBindings/` | Optional Python bindings |

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

**Prerequisites:** For RTX 50-series (`sm_120`), build **MFEM/dgtd with CUDA ≥ 13.x** (CUDA 12.8’s `cicc` can OOM on MFEM kernels; see [mfem#5363](https://github.com/mfem/mfem/issues/5363)). Set `-DCMAKE_CUDA_ARCHITECTURES` to your GPU (e.g. `120`, `89`, `86`) and match `CMakePresets.json`.

**HYPRE with CUDA** — stay on **HYPRE 2.31** (matches stock `mfem-geg`; no HYPRE-3 API patches). Build HYPRE with **CUDA 12.8** if needed (2.31 does not build against CUDA 13 CCCL); link that install into a dgtd/MFEM build that uses CUDA 13.3 for `sm_120`:
```sh
# Example: HYPRE 2.31 + CUDA 12.8, arch matched to the GPU
cmake -S hypre-2.31.0/src -B hypre-cuda-build \
      -DHYPRE_WITH_CUDA=ON \
      -DHYPRE_CUDA_SM=120 \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc \
      -DCMAKE_INSTALL_PREFIX=$HOME/workspace/hypre-cuda-install
cmake --build hypre-cuda-build -j $(nproc)
cmake --install hypre-cuda-build
export HYPRE_CUDA_DIR=$HOME/workspace/hypre-cuda-install
```

Then (with `METIS_DIR` set; CUDA 13.3 for the dgtd/MFEM compile via the preset):
```sh
export PATH=/usr/local/cuda-13.3/bin:$PATH
cmake --preset gnu-release-cuda
cmake --build --preset build-gnu-release-cuda --parallel
```

On `sm_120`, the top-level CMake compiles a small set of heavy MFEM CUDA TUs at `-O0` to keep `cicc` compile times reasonable.

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
./build/gnu-release-mpi/bin/opensemba_dgtd -i testData/maxwellInputs/1D_PEC/1D_PEC.json
```

CUDA builds (`gnu-release-cuda`, `gnu-release-cuda-sm120`, and the matching debug presets) default to `--device cuda`. Override with `--device cpu` or `--device omp` if you need the host. MPI-only binaries still default to `cpu`.

(Confirm binary name/path for your preset.)

Exports appear under `exports/` by kind, then run mode and case name:

- Simulation probes / stats / RCS dumps: `exports/SimulationData/<runmode>/<casename>/`
- ParaView exporter: `exports/ParaView/<runmode>/`
- Sparse operators: `exports/Operators/<casename>/`

## RCS post-processing (`opensemba_rcs`)

Offline far-field / RCS from an existing `rcssurface` export. The Maxwell case must already have been run with a `probes.rcssurface` probe so that

`exports/SimulationData/<runmode>/<casename>/RCSSurface/<probe_name>/rank*/surface_data.bin`

exists. The post-processor also reads the case JSON at `testData/maxwellInputs/<casename>/<casename>.json` (plane-wave data and probe names).

RCS input files live under `testData/rcsInputs/` (often `*.rcs.json`). Full parameter reference: **[docs/json-input-format.md](docs/json-input-format.md#offline-rcs-json-opensemba_rcs)**.

| Field | Required | Description |
|-------|----------|-------------|
| `runmode` | yes | Export folder segment (e.g. `cuda-1`, `single-core`, `mpi-8`) |
| `casename` | yes | Must match `testData/maxwellInputs/<casename>/` and the export tree |
| `frequencies` | yes | `{ start, end, steps }` in Hz (`steps` = linspace count) |
| `angles` | yes | `theta` / `phi` each `{ start, end, steps }` in radians |
| `max_time` | no | Keep snapshots with time ≤ this value |
| `every_n_steps` | no | Subsample while reading (default `1`; not equivalent to full history) |
| `ram_gate` | no | Hard RAM budget in **GiB**; omit → ~50% of `MemAvailable`. Over budget → streaming DFT |

Example:
```sh
mpiexec -n 1 ./build/gnu-release-mpi/bin/opensemba_rcs \
  -i testData/rcsInputs/3D_Nasa_Almond_G2_25cm_5GHz.rcs.json
```

Results are written under the same probe directory as `farfield/` and `rcs/`.

## Further documentation

| Topic | Link |
|-------|------|
| All docs | [docs/README.md](docs/README.md) |
| JSON input (solver + offline RCS) | [docs/json-input-format.md](docs/json-input-format.md) |
| MOR → ParaView | [docs/mor2paraview.md](docs/mor2paraview.md) |
| Box PML (mesh, JSON, absorption table) | [docs/designing-a-box-pml.md](docs/designing-a-box-pml.md) |

## Funding

- Spanish Ministry of Science and Innovation (MICIN/AEI) (Grant PID2022-137495OB-C31).
- European Union, HECATE project (HE-HORIZON-JU-Clean-Aviation-2022-01).
- European Union, FEDER 2020 (B-TIC-700-UGR20).
