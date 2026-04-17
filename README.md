[![ubuntu-gnu](https://github.com/OpenSEMBA/dgtd/actions/workflows/ubuntu-gnu.yml/badge.svg)](https://github.com/OpenSEMBA/dgtd/actions/workflows/ubuntu-gnu.yml)

# semba-dgtd
Maxwell's curl equations solver using discontinuous Galerkin methods

## Compiling

The project requires the following tools:
- CMake >= 3.25.2
- [vcpkg](https://github.com/microsoft/vcpkg) (pointed to by `$VCPKG_ROOT`)
- Ninja
- GCC

vcpkg will automatically install the following dependencies declared in [vcpkg.json](vcpkg.json):
- `eigen3`
- `gtest`
- `fftw3`
- `nlohmann-json`

### MPI builds

MPI builds additionally require METIS 5 and HYPRE to be compiled from sources and their install directories exported as environment variables before configuring.

**Build METIS:**
```sh
wget https://github.com/mfem/tpls/raw/gh-pages/metis-5.1.0.tar.gz
tar -zxvf metis-5.1.0.tar.gz
cd metis-5.1.0
make BUILDDIR=lib config
make BUILDDIR=lib
cp lib/libmetis/libmetis.a lib/
export METIS_DIR=$PWD
```

**Build HYPRE:**
```sh
wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.31.0.tar.gz
tar -zxvf v2.31.0.tar.gz
cd hypre-2.31.0/src
./configure --disable-fortran
make -j $(nproc)
export HYPRE_DIR=$PWD/hypre
```

### CUDA builds

CUDA builds require a separate HYPRE installation compiled with CUDA support. HYPRE's CMake build system is used for this (the autoconf `./configure` path does not support CUDA). The CUDA architecture must match the target GPU; the presets default to `sm_89` (Ada Lovelace / RTX 40-series).

**Build HYPRE with CUDA:**
```sh
wget https://github.com/hypre-space/hypre/archive/refs/tags/v2.31.0.tar.gz
tar -zxvf v2.31.0.tar.gz
cmake -S hypre-2.31.0/src -B hypre-cuda-build \
      -DHYPRE_WITH_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=89 \
      -DCMAKE_INSTALL_PREFIX=$HOME/hypre-cuda-install
cmake --build hypre-cuda-build -j $(nproc)
cmake --install hypre-cuda-build
export HYPRE_CUDA_DIR=$HOME/hypre-cuda-install
```

Then configure with a CUDA preset. `METIS_DIR` must also be set as described above:
```sh
cmake --preset gnu-release-cuda
cmake --build --preset build-gnu-release-cuda --parallel
```

### CMake presets

The project provides the following CMake presets:

| Preset | Description |
|--------|-------------|
| `gnu-debug-mpi` | Debug build with MPI and OpenMP |
| `gnu-release-mpi` | Release build with MPI and OpenMP |
| `gnu-debug-cuda` | Debug build with MPI, OpenMP and CUDA (gcc-12, sm_89) |
| `gnu-release-cuda` | Release build with MPI, OpenMP and CUDA (gcc-12, sm_89) |

Configure and build with a preset:
```sh
cmake --preset gnu-release-mpi
cmake --build --preset build-gnu-release-mpi --parallel
```

MFEM is included as a submodule (`external/mfem-geg`) and built automatically. To use an externally installed MFEM instead, set `-DSEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY=OFF`.

> **Warning:** The required MFEM installation must be the [OpenSEMBA/mfem-geg](https://github.com/OpenSEMBA/mfem-geg) fork. This repository contains modifications that are not present in upstream MFEM; the project will not compile against a standard MFEM installation.

### Defining a JSON file

Cases can be defined by parsing a JSON file with the problem information. See a complete [example](https://github.com/OpenSEMBA/dgtd/blob/main/testData/maxwellInputs/1D_PEC/1D_PEC.json) of a valid JSON file.
An OpenSEMBA/dgtd JSON file must have the following structure; **bold** entries are **required**:

- solver_options:
	- Object. User can customise solver settings. If undefined, all defaults apply.
		- evolution_operator:
			- String. Selects the DG evolution operator. Can be `"maxwell"`, `"global"`, or `"hesthaven"`. If undefined, defaults to `"global"`.
		- upwind_alpha:
			- Double. Upwind flux blending factor. `0.0` = fully centered, `1.0` = fully upwind. If undefined, defaults to `1.0`.
		- final_time:
			- Double. Total simulation duration in natural units (1 meter/c). If undefined, defaults to `2.0`.
		- time_step:
			- Double. Fixed time step in natural units. Must be defined for 2D and 3D problems. Overrides `cfl` for 1D if both are defined. If undefined or `0.0` in 1D, an automatic step is computed via CFL.
		- cfl:
			- Double. Courant–Friedrichs–Lewy condition used to compute the automatic time step in 1D. Not available for 2D or 3D. Ignored if `time_step` is also defined. Must be in the range (0.0, 1.0].
		- order:
			- Integer. Polynomial order of the finite element basis. If undefined, defaults to `3`.
		- spectral:
			- Boolean. Use a spectral evolution operator that assembles the full E/H system matrix and derives the time step from its eigenvalues. High computational cost; does not support all features. If undefined, defaults to `false`.
		- export_operator:
			- Boolean. Write the assembled evolution operator matrix to disk for inspection. If undefined, defaults to `false`.
		- basis_type:
			- String. Finite element basis type passed to MFEM.
		- ode_type:
			- String. ODE time-integration method passed to MFEM.

- **model**:
	- Object. Contains geometry, material, and boundary information.
		- **filename**:
			- String. Mesh filename (`.msh` or `.mesh`). The file must reside in the same directory as the JSON file.
		- refinement:
			- Integer. Number of uniform refinement levels to apply to the loaded mesh. Optional.
		- **materials**:
			- Array. At least one entry is required. Each entry defines the electromagnetic properties of one or more mesh domains:
				- **tags**:
					- Array of integers. Mesh volume/surface/segment attribute IDs that share these material properties. (Volumes in 3D, Surfaces in 2D, Segments in 1D.)
				- relative_permittivity:
					- Double. Relative permittivity $\varepsilon_r$. Defaults to `1.0`.
				- relative_permeability:
					- Double. Relative permeability $\mu_r$. Defaults to `1.0`.
				- bulk_conductivity:
					- Double. Bulk electrical conductivity in S/m. Defaults to `0.0`. Internally scaled by the free-space impedance.
		- **boundaries**:
			- Array. At least one entry is required. Each entry defines a boundary or interface condition:
				- **tags**:
					- Array of integers. Mesh boundary attribute IDs. (Surfaces in 3D, Segments in 2D, Points in 1D.)
				- **type**:
					- String. Boundary condition class. Can be `"PEC"`, `"PMC"`, `"SMA"`, or `"SGBC"`.
				- *(For `type: "SGBC"` only)* material:
					- Object. Defines a single-layer surface general boundary condition. Mutually exclusive with `layers`.
						- relative_permittivity: Double. Defaults to `1.0`.
						- relative_permeability: Double. Defaults to `1.0`.
						- **bulk_conductivity**: Double. Conductivity of the layer material in S/m.
						- **material_width**: Double. Physical thickness of the layer in meters.
						- num_of_segments: Integer. Number of sub-elements in the layer mesh. Auto-computed from material properties and source frequency if omitted.
						- order: Integer. Polynomial order for the layer sub-solver. Auto-computed if omitted.
				- *(For `type: "SGBC"` only)* layers:
					- Array. Defines a multi-layer SGBC stack. Mutually exclusive with `material`. Each element has the same fields as `material` above.
				- *(For `type: "SGBC"` only)* sgbc_boundaries:
					- Object. Boundary conditions applied to the outer and inner faces of the SGBC sub-domain. Can be placed at the boundary level or inside `material` (the latter is kept for backward compatibility).
						- left: String. Condition on the outer (field-side) face. Can be `"PEC"`, `"PMC"`, or `"SMA"`.
						- right: String. Condition on the inner face. Can be `"PEC"`, `"PMC"`, or `"SMA"`.

- probes:
	- Object. Controls data extraction. If undefined, no data is recorded.
		- exporter:
			- Object. Enables ParaView (VisIt) field export.
				- name: String. Output dataset name. Defaults to the mesh filename stem.
				- steps: Integer. Export every N time steps. Mutually exclusive with `saves`.
				- saves: Integer. Total number of exports over the whole simulation. The step interval is computed automatically. Mutually exclusive with `steps`.
		- point:
			- Array. Each entry records all E and H field components at a single point every interval.
				- **position**: Array of doubles. Spatial coordinates. Must match the mesh dimension (e.g. `[x, y]` for 2D). ***Warning:*** If the point lies outside the mesh the simulation will crash.
				- steps / saves: See exporter above.
		- field:
			- Array. Each entry records a single scalar field component at a point.
				- **field_type**: String. Can be `"electric"` or `"magnetic"`.
				- **polarization**: String. Component to record. Can be `"X"`, `"Y"`, or `"Z"`.
				- **position**: Array of doubles. Spatial coordinates. Same constraints as for `point`.
				- steps / saves: See exporter above.
		- farfield:
			- Array. Each entry accumulates near-to-far-field data for 3D RCS post-processing.
				- **tags**: Array of integers. Mesh surface tags that form the near-field Huygens surface.
				- name: String. Output name.
				- export_path: String. Directory for the output files.
				- steps / saves: See exporter above.
		- domain_snapshot:
			- Object. Periodic full-domain field snapshot (alternative to the incremental exporter).
				- name: String. Output name. Defaults to the mesh filename stem.
				- steps / saves: See exporter above.
		- rcssurface:
			- Array. Each entry performs in-situ 2D RCS surface integration.
				- **tags**: Array of integers. Mesh surface tags for the integration surface.
				- name: String. Output name.
				- steps / saves: See exporter above.

- **sources**:
	- Array. Defines the electromagnetic excitation. At least one source is required.
		- **type**: String. Can be `"initial"`, `"planewave"`, or `"dipole"`.

		*For `type: "initial"` — volumetric initial condition:*
		- **field_type**: String. Can be `"electric"` or `"magnetic"`.
		- **polarization**: Array of 3 doubles. Polarization direction vector.
		- **magnitude**: Object.
			- **type**: String. Can be `"gaussian"`, `"resonant"`, `"besselj6_2D"`, or `"besselj6_3D"`.
			- *(For `type: "gaussian"`)* **spread**: Double. Standard deviation $\sigma$ of the Gaussian pulse.
			- *(For `type: "resonant"`)* **modes**: Array of integers. Number of standing waves along each spatial axis.
		- *(Required for `magnitude.type: "gaussian"`)* **center**: Array of n doubles. Spatial position of the Gaussian centroid.
		- *(Required for `magnitude.type: "gaussian"`)* **dimension**: Integer. Number of active spatial dimensions in the Gaussian exponent.

		*For `type: "planewave"` — total-field/scattered-field (TFSF) plane wave:*
		- **tags**: Array of integers. Mesh boundary tags that form the TFSF interface surface.
		- **polarization**: Array of 3 doubles. E-field polarization direction.
		- **propagation**: Array of 3 doubles. Wave propagation direction vector.
		- **magnitude**: Object.
			- **spread**: Double. Standard deviation $\sigma$ of the Gaussian envelope.
			- **mean**: Array of doubles. Position of the Gaussian center projected onto the propagation axis at $t = 0$. Typically set behind the TFSF box. Provide as a 1-, 2-, or 3-component vector matching the mesh dimension.
			- frequency: Double (Hz). Carrier frequency. If defined, wraps the Gaussian in a sinusoidal modulation (modulated Gaussian). If omitted, a broadband Gaussian pulse is used.

		*For `type: "dipole"` — derivative-Gaussian dipole source within a TFSF box:*
		- **tags**: Array of integers. Mesh boundary tags that form the TFSF interface surface.
		- **magnitude**: Object.
			- **length**: Double. Dipole length.
			- **spread**: Double. Gaussian spread parameter.
			- **mean**: Double. Gaussian center position along the dipole axis.

## Funding

- Spanish Ministry of Science and Innovation (MICIN/AEI) (Grant Number: PID2022-137495OB-C31).
- European Union, HECATE project (HE-HORIZON-JU-Clean-Aviation-2022-01).
- European Union, FEDER 2020 (B-TIC-700-UGR20).
