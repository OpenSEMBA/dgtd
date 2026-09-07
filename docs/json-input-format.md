# JSON input format

Simulation cases are defined by JSON files parsed at runtime. Each case must follow the repository naming layout: the case folder name, the `.json` filename, the `.msh` filename, and `model.filename` must all share the same `<case_name>` base, and `model.filename` must not include a path.

Example: [testData/maxwellInputs/1D_PEC/1D_PEC.json](../testData/maxwellInputs/1D_PEC/1D_PEC.json)

Legend: **[REQUIRED]** = must be present; **[OPTIONAL]** = has a default value (shown).

## solver_options

Object. User can customise solver settings. If undefined, all defaults apply.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `evolution_operator` | string | `"global"` | DG evolution operator: `"maxwell"`, `"global"`, or `"hesthaven"`. **New features should use `"global"` only.** |
| `upwind_alpha` | double | `1.0` | Upwind flux blending: `0.0` = centered, `1.0` = upwind. |
| `final_time` | double | `2.0` | Simulation duration in natural units (1 m/c). |
| `time_step` | double | `0.0` | Fixed time step in natural units. Required for 2D/3D. In 1D, `0.0` triggers automatic CFL-based step. |
| `cfl` | double | `1.0` | CFL for automatic 1D time step. Range (0, 1]. Ignored if `time_step` is set. Not used in 2D/3D auto step. |
| `sgbc_cfl` | double | `0.5` | Crossing-time CFL for SGBC sub-step recommendation: `δt_rec = crossing_time × sgbc_cfl × opacity_relax` (then SI→natural). Must be `> 0`. Omit for historical default (`0.5`). Does not change SGBC layer geometry; only `nsteps = ceil(Δt / δt_rec)`. |
| `order` | integer | `2` | Polynomial order of the FE basis. |
| `spectral` | boolean | `false` | Spectral evolution operator (full matrix eigenvalue step). High cost; limited feature support. |
| `export_operator` | boolean | `false` | Write assembled evolution operator to disk. |
| `basis_type` | integer | `1` | MFEM basis: `0` GaussLegendre, `1` GaussLobatto, `2` Bernstein, `3` OpenUniform, `4` CloseUniform, `5` OpenHalfUniform. |
| `ode_type` | integer | `0` | Time integrator: `0` RK4, `1` BackwardEuler, `2` Trapezoidal, `3` ImplicitMidpoint, `4` SDIRK33, `5` SDIRK23, `6` SDIRK34. |

### evolution_operator: hesthaven (fast explicit path)

`"hesthaven"` uses element-local dense operators (matrix-free on straight meshes). Use for explicit RK4 cases **without** SGBC, volumetric PML, bulk conductivity, or implicit `ode_type`.

| Capability | hesthaven | global |
|------------|-----------|--------|
| MPI (`mpirun`) | Yes — shared-face ghost exchange and neighbor connectivity in `Mult()` | Yes |
| CUDA (`--device cuda`) | Yes — element kernels when build has `SEMBA_DGTD_ENABLE_CUDA` | Yes |
| SGBC / PML / conductivity | No | Yes |
| Implicit `ode_type` | No | Yes |
| Centered SMA (`upwind_alpha: 0`) | Blocked in driver | Yes |

Run example: `mpirun -np 4 ./opensemba_dgtd case.json --device cuda` with `"evolution_operator": "hesthaven"`.

## model

Object. Geometry, materials, and boundaries.

### filename [REQUIRED]

String. Mesh file (`.msh` or `.mesh`) in the same directory as the JSON file.

### refinement [OPTIONAL]

Integer. Uniform mesh refinement levels.

### materials [REQUIRED]

Array. At least one entry. Each entry assigns electromagnetic properties to mesh domains (volumes in 3D, surfaces in 2D, segments in 1D).

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tags` | int[] | — | Mesh attribute IDs sharing these properties. |
| `type` | string | (legacy) | `"vacuum"` or `"PML"` (volumetric ADE-PML region — see [pml/03-json-schema-and-mesh.md](./pml/03-json-schema-and-mesh.md); CFS paused). If omitted, legacy eps/mu/sigma fields apply. |
| `relative_permittivity` | double | `1.0` | ε_r (legacy / non-PML). |
| `relative_permeability` | double | `1.0` | μ_r (legacy / non-PML). |
| `bulk_conductivity` | double | `0.0` | Conductivity in S/m; scaled internally by free-space impedance. Not for PML tags. |

### boundaries [REQUIRED]

Array. At least one entry per boundary or interface condition.

| Field | Type | Description |
|-------|------|-------------|
| `tags` | int[] | Boundary attribute IDs. |
| `type` | string | `"PEC"`, `"PMC"`, `"SMA"`, or `"SGBC"`. |

#### type: SGBC

Surface general boundary condition (multi-layer absorber / material slab sub-solver).

| Field | Description |
|-------|-------------|
| `exporter_probe` | If `true` and top-level `probes.exporter` exists, export internal SGBC fields as `InsideSGBC_tag<first_tag>`. |
| `material` | Single layer (mutually exclusive with `layers`). |
| `layers` | Multi-layer stack (mutually exclusive with `material`). |

Layer / material object fields:

| Field | Default | Description |
|-------|---------|-------------|
| `relative_permittivity` | `1.0` | |
| `relative_permeability` | `1.0` | |
| `bulk_conductivity` | — | S/m |
| `material_width` | — | Layer thickness (m) |
| `num_of_segments` | auto | Sub-elements in layer mesh |
| `order` | auto | Polynomial order for layer sub-solver |

`sgbc_boundaries` (SGBC only): inner/outer face conditions.

| Field | Values |
|-------|--------|
| `left` | `"PEC"`, `"PMC"`, `"SMA"` (field-side) |
| `right` | `"PEC"`, `"PMC"`, `"SMA"` (inner face) |

## probes [OPTIONAL]

If omitted, no probe output.

### exporter

ParaView (VisIt) field export.

| Field | Description |
|-------|-------------|
| `name` | Dataset name (default: mesh basename) |
| `steps` | Export every N steps (exclusive with `saves`) |
| `saves` | Total exports over run (interval computed automatically) |

### point

Array. All E/H components at a point.

| Field | Description |
|-------|-------------|
| `position` | Coordinates matching mesh dimension |
| `steps` / `saves` | As above |

**Warning:** point outside mesh crashes the simulation.

### field

Array. Single scalar component at a point.

| Field | Values |
|-------|--------|
| `field_type` | `"electric"`, `"magnetic"` |
| `polarization` | `"X"`, `"Y"`, `"Z"` |
| `position` | Coordinates |
| `steps` / `saves` | As above |

### farfield

Near-to-far-field surface export under `Exports/<run-mode>/<case>/NearToFarFieldProbes/<name>/`.

| Field | Description |
|-------|-------------|
| `tags` | Boundary tags for NFS surface |
| `name` | Default `"NearFieldProbe"` |
| `steps` / `saves` | As above |

### domain_snapshot

Periodic full-domain snapshot (alternative to incremental exporter).

### rcssurface

Surface E/H snapshots for offline RCS (`surface_data.bin` with geometry header + time blocks).

| Field | Description |
|-------|-------------|
| `tags` | Surface integration tags |
| `name` | Default `"RCSSurfaceProbe"` |
| `steps` / `saves` | As above. Prefer matching temporal density across cases (e.g. with `time_step` 0.0001 use `steps: 5` to match a 0.0005 / `steps: 1` export). |

Offline `opensemba_rcs` input JSON may also set `every_n_steps` (integer ≥ 1, default 1) to subsample an existing dense `surface_data.bin` **while reading** (skips payloads; does not load them into RAM).

### mor_state

Full DG state vector snapshots compatible with exported `{name}_global.csr` for offline `y = A x`.

| Field | Description |
|-------|-------------|
| `record_time_start` | Recording start time |
| `record_time_final` | Recording end time |
| `saves` | Number of uniform snapshots |
| `name` | Default `"MORState"` |

Snapshot file format (`x_0`, `x_1`, …): line 1 = time, line 2 = size, then DOF values in `[Ex | Ey | Ez | Hx | Hy | Hz]` order.

## sources [REQUIRED]

Array. At least one source; all entries superimpose.

| `type` | Description |
|--------|-------------|
| `"initial"` | Volumetric initial condition |
| `"planewave"` | TFSF plane wave |
| `"dipole"` | TFSF dipole |

### type: initial

| Field | Description |
|-------|-------------|
| `field_type` | `"electric"` or `"magnetic"` |
| `polarization` | 3-vector |
| `magnitude.type` | `"gaussian"`, `"resonant"`, `"besselj6_2D"`, `"besselj6_3D"` |
| `magnitude.spread` | Gaussian σ (for gaussian) |
| `magnitude.modes` | Standing-wave modes (for resonant) |
| `center` | Gaussian centroid (required for gaussian) |
| `dimension` | Active dimensions in Gaussian (required for gaussian) |

### type: planewave

| Field | Description |
|-------|-------------|
| `tags` | TFSF interface boundary tags |
| `polarization` | E polarization |
| `propagation` | Propagation direction |
| `magnitude.spread` | Gaussian envelope σ |
| `magnitude.mean` | Optional pulse center on propagation axis |
| `magnitude.frequency` | Optional carrier (Hz) for modulated Gaussian |

### type: dipole

| Field | Description |
|-------|-------------|
| `tags` | TFSF interface tags |
| `magnitude.length` | Dipole length |
| `magnitude.spread` | Gaussian spread |
| `magnitude.mean` | Optional center along dipole axis |
