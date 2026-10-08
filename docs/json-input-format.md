# JSON input format

Simulation cases are defined by JSON files parsed at runtime. Each case must follow the repository naming layout: the case folder name, the `.json` filename, the `.msh` filename, and `model.filename` must all share the same `<case_name>` base, and `model.filename` must not include a path.

Example: [testData/maxwellInputs/1D_PEC/1D_PEC.json](../testData/maxwellInputs/1D_PEC/1D_PEC.json)

Legend: **[REQUIRED]** = must be present; **[OPTIONAL]** = has a default value (shown).

## solver_options

Object. User can customise solver settings. If undefined, all defaults apply.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `evolution_operator` | string | `"global"` | `"global"` (default; SGBC, volumetric PML, Debye, Lorentz, conductivity). `"hesthaven"` is a limited explicit operator without those features. `"maxwell"` is deprecated. |
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
| `checkpoint_percent` | double | `20` | `0` disables checkpoints. Otherwise, after the completed step that crosses each multiple of this percent of `final_time`, write the evolved state. The default `20` saves near 20%, 40%, 60%, 80%, and 100%. |

### Checkpoints

`checkpoint_percent` writes under `exports/SimulationData/<run-mode>/<case>/Checkpoints/`. Each save stores the packed field state (E, H, and any PML, Debye, or Lorentz auxiliaries), SGBC face history, the time and time step, probe cursors, the element partition, and the solver run time up to that save. Operators are rebuilt on the next launch. A finished run deletes that directory. Probe files, ParaView output, and `SimulationStats` stay.

`Simulation Run Time` in `statistics_rank*.dat` is the time spent inside the time loop and covered by a checkpoint, plus the time after the last resume. Minutes after the last save and before a crash are not included. A run that reached its 40% save in 4 minutes, died at 50%, and then took 6 minutes to finish after `--restart` reports 10 minutes.

Resume with the same JSON, the same mesh, and the same MPI rank count and device, from the same working directory:

```sh
mpirun -np N ./opensemba_dgtd -i case.json --restart
```

`solver_options.checkpoint_percent` may differ on resume. Any other JSON change, or a different mesh, is rejected. A checkpoint written with another `-np` or device is rejected as well, including when `--restart` is omitted, so that launch does not start at t = 0. Delete the `Checkpoints` directory to start again at t = 0. `SIGINT` or `SIGTERM` delivered to the `opensemba_dgtd` processes finishes the current step, writes one checkpoint, and exits. `mpirun` handles Ctrl-C itself and can kill those processes before that write; the last complete checkpoint is what a later `--restart` loads. A failed save does not stop the run. The previous complete checkpoint is kept, and the next attempt waits at least a minute.

### evolution_operator: hesthaven

`"hesthaven"` is element-local dense operators (matrix-free on straight meshes). It does **not** run SGBC, volumetric PML, Debye, Lorentz, bulk conductivity, or implicit `ode_type`. Use `"global"` (or omit `evolution_operator`) for those features.

| Capability | hesthaven | global |
|------------|-----------|--------|
| MPI (`mpirun`) | Yes — shared-face ghost exchange and neighbor connectivity in `Mult()` | Yes |
| CUDA | Yes — CUDA builds of `opensemba_dgtd` default to `--device cuda`; pass `--device cpu` (or `omp`) to stay on the host | Yes |
| SGBC / PML / Debye / Lorentz / conductivity | No | Yes |
| Implicit `ode_type` | No | Yes |
| Centered SMA (`upwind_alpha: 0`) | Blocked in driver | Yes |

Run example: `mpirun -np 4 ./opensemba_dgtd -i case.json` with `"evolution_operator": "hesthaven"` (CUDA binaries default to the GPU).

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
| `type` | string | (legacy) | `"vacuum"` or `"PML"`. PML fields, mesh tagging, and the absorption reference: [designing-a-box-pml.md](./designing-a-box-pml.md). `"R"`/`"r"`, `sigma_max`, `stretch_mode: "radial"`, and `radial_center` are rejected at parse. |
| `relative_permittivity` | double | `1.0` | ε_r (legacy / non-PML). |
| `relative_permeability` | double | `1.0` | μ_r (legacy / non-PML). |
| `bulk_conductivity` | double | `0.0` | Conductivity in S/m; scaled internally by free-space impedance. Not for PML tags. |
| `debye` | object | — | Single-pole electric Debye on an untyped material. See below. |
| `lorentz` | object | — | Single-pole electric Lorentz on an untyped material. See below. |

#### debye

Optional object. Legal only when the material entry has no `type`. The electric mass uses `eps_inf`. `tau` is seconds; the solver stores $\tau_{\mathrm{SI}}\,c_{\mathrm{SI}}$ because time is in light-meters.

| Field | Type | Description |
|-------|------|-------------|
| `eps_inf` | double | Instantaneous relative permittivity. Must be $\ge 1$. |
| `eps_s` | double | Static relative permittivity. Must be greater than `eps_inf`. |
| `tau` | double | Relaxation time in seconds. Must be $> 0$. |

`relative_permittivity` is rejected together with `debye`. `relative_permeability` and `bulk_conductivity` keep their usual defaults. `bulk_conductivity` stays an independent Ohm term beside the pole.

`debye` is rejected on `type: "vacuum"`, `type: "PML"`, and on an SGBC layer. Initialization also aborts if a Debye tag is a PML tag. A tag cannot be listed as both Debye and another material. `evolution_operator` must be `"global"`. `spectral: true` is rejected. Implicit ODE types abort when the state includes the Debye polarization.

#### lorentz

Optional object. Legal only when the material entry has no `type`. The electric mass uses `eps_inf`. Give the three rates in **one** style, never mixed. The solver stores $\omega/c_{\mathrm{SI}}$ because time is in light-meters. $f_1 = 0$ (or $\omega_1 = 0$) is a cold plasma on this same pole.

The relative permittivity is

$$
\varepsilon_r(\omega)=\varepsilon_\infty+\frac{\omega_p^2}{\omega_1^2-\omega^2+j\,2\gamma\omega}.
$$

**Hz style (preferred for new cases):** `eps_inf`, `f_p`, `f_1`, `gamma`. All three rates are ordinary frequencies in Hz. Conversion is $\omega=2\pi f$, including damping: JSON `gamma: 2.5e7` is $25\,\mathrm{MHz}$, stored as $\gamma=2\pi\times 2.5\times 10^7\,\mathrm{rad/s}$.

**rad/s style:** `eps_inf`, `omega_p`, `omega_1`, `gamma`. All three rates are angular frequencies in rad/s. Do not mix with `f_p` / `f_1`.

| Field | Type | Description |
|-------|------|-------------|
| `eps_inf` | double | Instantaneous relative permittivity. Must be $\ge 1$. |
| `f_p` | double | Plasma frequency in Hz. Must be $> 0$. |
| `f_1` | double | Resonance frequency in Hz. Must be $\ge 0$. |
| `omega_p` | double | Plasma frequency in rad/s. Must be $> 0$. |
| `omega_1` | double | Resonance frequency in rad/s. Must be $\ge 0$. |
| `gamma` | double | Damping. Hz with `f_*`, or rad/s with `omega_*`. Must be $\ge 0$. |

Example (Hz), from [3D_RCS_Lorentz_G2](../testData/maxwellInputs/3D_RCS_Lorentz_G2/3D_RCS_Lorentz_G2.json):

```json
"lorentz": {
  "eps_inf": 2.0,
  "f_p": 2.5e8,
  "f_1": 2.0e8,
  "gamma": 2.5e7
}
```

That is the same pole as `omega_p = 2π f_p`, `omega_1 = 2π f_1`, `gamma = 2π × 2.5e7` rad/s. Older 1D cases still use the rad/s keys; both styles are accepted.

`relative_permittivity` is rejected together with `lorentz`. `bulk_conductivity` stays an independent Ohm term. A tag cannot carry both `debye` and `lorentz`. `lorentz` is rejected on vacuum, PML, and SGBC. Initialization also aborts if a Lorentz tag is a PML tag. `evolution_operator` must be `"global"`. `spectral: true` is rejected. Implicit ODE types abort when the state includes Lorentz auxiliaries.

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

ParaView (VisIt) field export under `exports/ParaView/<run-mode>/`.

| Field | Description |
|-------|-------------|
| `name` | Dataset name (default: mesh basename) |
| `save_every` | **Required cadence.** Export interval in solver time. Snapshots at `t = 0`, `save_every`, `2·save_every`, … and **always** at `final_time`. Example: `final_time: 20` with `save_every: 0.5` → 41 frames (`0, 0.5, …, 20`). Must be `> 0`. |

Legacy keys on `exporter`:
- `saves` (total frame count) — **rejected**; convert with `save_every = final_time / (saves - 1)` when `saves > 1`.
- `steps` (every N time steps) — still parsed for backward compatibility, but prefer `save_every = steps * time_step`. Do not combine with `save_every`.

### point

Array. All E/H components at a point.

| Field | Description |
|-------|-------------|
| `position` | Coordinates matching mesh dimension |
| `steps` | Record every N time steps (exclusive with `saves`) |
| `saves` | Total samples over the run (step interval computed from `final_time` / `time_step`) |

**Warning:** point outside mesh crashes the simulation.

### field

Array. Single scalar component at a point.

| Field | Values |
|-------|--------|
| `field_type` | `"electric"`, `"magnetic"` |
| `polarization` | `"X"`, `"Y"`, `"Z"` |
| `position` | Coordinates |
| `steps` / `saves` | Same meaning as `point` (count / step interval; not `save_every`) |

### farfield

Near-to-far-field surface export under `exports/SimulationData/<run-mode>/<case>/NearToFarFieldProbes/<name>/`.

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

Offline frequency/angle sweeps use a **separate** JSON for `opensemba_rcs` — see [Offline RCS JSON](#offline-rcs-json-opensemba_rcs). Do not confuse offline `every_n_steps` with case-JSON `probes.rcssurface.steps` (export cadence during `opensemba_dgtd`).

### mor_state

Full DG state vector snapshots compatible with exported `{name}_global.csr` for offline `y = A x`. Written under `exports/SimulationData/<run-mode>/<case>/MORStateProbes/<name>/`.

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
| `"delta_gap"` | Impressed tangential E on an interior face (2D curve or 3D surface) |
| `"coaxial_port"` | TFSF coaxial TEM mode on an interior annular face |

### Gaussian width (`spread` / `f_1e`)

Time-domain Gaussians are $\exp(-(t-t_0)^2/(2\sigma^2))$ with $\sigma$ in light-metres. Set the width with **either**:

- `spread`: $\sigma$ directly, or
- `f_1e`: frequency in Hz where incident **power** is $1/e$ of DC. Then $\sigma = c/(2\pi f_{1e})$.

Do not mix units: $f_{1e}=3\times 10^8$ is 300 MHz, which is $\sigma\approx 0.16$, not `spread: 0.15` written as a frequency. If both keys are present, `f_1e` is used and a warning is printed. At solver start that 1/e frequency is compared to mean element size and FE order; if the pulse is likely under-resolved a warning is printed and the run continues.

On `planewave`, `dipole`, and `initial`, the keys live under `magnitude`. On `delta_gap` and `coaxial_port`, they are on the source object.

### type: initial

| Field | Description |
|-------|-------------|
| `field_type` | `"electric"` or `"magnetic"` |
| `polarization` | 3-vector |
| `magnitude.type` | `"gaussian"`, `"resonant"`, `"besselj6_2D"`, `"besselj6_3D"` |
| `magnitude.spread` | Gaussian $\sigma$ in light-metres |
| `magnitude.f_1e` | Alternative to `spread`: Hz at 1/e incident power |
| `magnitude.modes` | Standing-wave modes (for resonant) |
| `center` | Gaussian centroid (required for gaussian) |
| `dimension` | Active dimensions in Gaussian (required for gaussian) |

### type: planewave

| Field | Description |
|-------|-------------|
| `tags` | TFSF interface boundary tags |
| `polarization` | E polarization |
| `propagation` | Propagation direction |
| `magnitude.spread` | Gaussian envelope $\sigma$ in light-metres |
| `magnitude.f_1e` | Alternative to `spread`: Hz at 1/e incident power |
| `magnitude.mean` | Optional pulse center on propagation axis |
| `magnitude.frequency` | Optional carrier (Hz) for modulated Gaussian |

Example:

```json
"magnitude": { "f_1e": 3.0e8 }
```

### type: dipole

| Field | Description |
|-------|-------------|
| `tags` | TFSF interface tags |
| `magnitude.length` | Dipole length (Hertzian geometric scale) |
| `magnitude.spread` | Gaussian spread σ in light-metres. Alternative: `f_1e` (Hz, 1/e incident power) |
| `magnitude.amplitude_peak` | Desired max equatorial \|E\| (\|E_θ\|) at `peak_radius` (default `1.0`). The analytic formula is scaled so that peak equals this value. |
| `magnitude.peak_radius` | Radius used to define `amplitude_peak` (default: min radius on the dipole TFSF tags if available, else `1.0`). |
| `magnitude.mean` | Optional center along dipole axis |

### type: delta_gap

One entry, on a 2D or 3D mesh. The tagged face must be interior, with vacuum on both sides. The source is a sibling of TFSF: its own face matrix, not a `planewave` or `dipole`.

Peak field is `magnitude` (default `1.0`). The pulse width is `spread` (Gaussian $\sigma$ in light-metres) or `f_1e` (Hz where incident power is $1/e$ of DC). The pulse is $e^{-(t-t_0)^2/(2\sigma^2)}$, and its center is five widths after $t = 0$. If both `spread` and `f_1e` are set, `f_1e` wins and a warning is printed.

If both `spread` and `f_1e` are omitted, $\sigma$ is chosen from the gap length $L$ so the spectrum is `db_cut` (default `-20`) down at $\mathrm{magnitude} \times L / 10$. That product is not a frequency you set.

Optional `signal` is `"gaussian"` (default) or `"gaussian_derivative"`. The derivative is what a flux must use if that pulse is the field: the time integral of the derivative is the pulse, so the slot charges and then empties. `spread` and `magnitude` describe that pulse in both cases. `magnitude` remains the peak $|E|$.

`polarization` is a required nonzero 3-vector, the direction of the electric field. It is normalized, and it must lie in the gap face. `magnitude` remains the peak $|E|$. For the vertical gap whose normal is $+x$, the field across the slot is $[0, 1, 0]$.

`magnitude` is the peak $|E^{\mathrm{inc}}|$ (default `1.0`), the zero-thickness stand-in for a gap voltage over a gap width. The magnetic incident field is zero. Both sides of the face carry the same electric trace, $+\tfrac{1}{2}$ times that polarization. The source keeps each side's own face block, so those equal traces do not cancel.

In 3D the tagged face is a rectangle. $h$ is the span of its vertices along `polarization`, the plate separation. $w$ is the span in the face perpendicular to that direction, the plate width. With the wave impedance equal to 1, the wide-plate line impedance is $Z=h/w$ and the capacitance per unit length along the face normal is $C'=w/h$. The run prints both, together with $w/h$. They are the parallel-plate consequences of that rectangle. They do not scale the source, and $C'$ is not a lumped capacitor: the face has no thickness along the normal. The approximation drops fringing, so it wants $w\gg h$.

| Field | Description |
|-------|-------------|
| `tags` | Interior face tags. A curve in 2D, a surface in 3D |
| `polarization` | Required 3-vector. Direction of $E$. Must be tangent to the gap |
| `magnitude` | Peak \|E\|. Default `1.0` |
| `spread` | Gaussian width $\sigma$ in light-metres. Must be `> 0`. Default from $L$ and `db_cut` |
| `f_1e` | Alternative to `spread`: Hz at 1/e incident power |
| `db_cut` | Used only when `spread` and `f_1e` are omitted. Spectrum level in dB at $\mathrm{magnitude} \times L / 10$. Default `-20`. Must be negative |
| `signal` | `"gaussian"` (default) or `"gaussian_derivative"` |

### type: coaxial_port

One entry, on a 3D mesh. An interior load face is the TFSF interface: the wave is launched into the total-field volume, and the scattered-field volume carries no incident field. A boundary load face has no second volume. It is total field only, the tag must also be SMA, and the axis points from that face into the single volume. The face then applies the SMA flux to the solution minus the TEM pair. Do not combine it with `planewave` or `dipole`.

`tags.outer` and `tags.live` are the PEC conductors used to measure the shared center and the radii. Cylindrical tags are kept. End caps, whose radius is not constant, are ignored. `tags.load` is the annular face. It is not a boundary condition.

The scattered-field volume is the load neighbor that meets an SMA boundary when exactly one side does. An interior SMA face meets both volumes on that face. When both sides meet an SMA boundary, the scattered-field volume is Elem1 of the load face and the total-field volume is Elem2. Every load face must share that order. When neither side meets an SMA boundary, it is the load neighbor that reaches a PML volume without crossing the load. If both sides reach a PML, the scattered-field side is the one whose PML interface is closer to the load. The other load neighbor is the total-field volume. The existing TFSF face operator and its Elem1/Elem2 convention are used as they are.

The incident field is the circular TEM pair, with magnitude the voltage $V_0$:

$$
E_\rho=\frac{V_0}{\rho\ln(b/a)},\qquad \mathbf{H}=\hat{\mathbf{s}}\times\mathbf{E}.
$$

$\hat{\mathbf{s}}$ points from the scattered-field volume into the total-field volume. In solver units the cable impedance is $\ln(b/a)/(2\pi)$.

`magnitude` defaults to `1.0` and is $V_0$, not a uniform peak $|E|$. Pulse width is `spread` (default `1.0`, light-metres) or `f_1e` (Hz, 1/e incident power). If both are set, `f_1e` wins and a warning is printed. The pulse center is five widths after $t=0$, the same delay as `delta_gap`. `signal` is `"gaussian"` (default) or `"gaussian_derivative"`.

| Field | Description |
|-------|-------------|
| `tags.outer` | Outer-conductor surface tags |
| `tags.live` | Inner-conductor surface tags |
| `tags.load` | Interior annular face. Becomes the TFSF marker |
| `magnitude` | Voltage $V_0$. Default `1.0` |
| `spread` | Gaussian width $\sigma$ in light-metres. Default `1.0` |
| `f_1e` | Alternative to `spread`: Hz at 1/e incident power |
| `signal` | `"gaussian"` (default) or `"gaussian_derivative"` |

---

## Offline RCS JSON (`opensemba_rcs`)

Separate input file for the offline RCS / far-field post-processor (`opensemba_rcs`). It does **not** replace the Maxwell case JSON; it points at an existing export and at `testData/maxwellInputs/<casename>/<casename>.json`.

**Prerequisites**

1. Case JSON includes `probes.rcssurface` (see [rcssurface](#rcssurface)).
2. `opensemba_dgtd` has been run so that  
   `exports/SimulationData/<runmode>/<casename>/RCSSurface/<probe_name>/rank*/surface_data.bin`  
   (and `mesh`) exist.
3. Launch from the repository root (paths are relative: `./exports/...`, `./testData/maxwellInputs/...`).

**Process**

1. Resolve `casename` → case JSON; collect each `rcssurface` probe `name`.
2. For each probe, read `exports/SimulationData/<runmode>/<casename>/RCSSurface/<name>/`.
3. Estimate load-all peak RAM vs budget (`ram_gate` or ~50% of Linux `MemAvailable`).
4. If the dump fits: load all kept snapshots and DFT. If not: stream snapshots into the frequency-domain accumulator (same kept samples ⇒ near-equivalent RCS).
5. NTFF / RCS over the requested frequency and angle grids; write `farfield/` and `rcs/` under the probe directory.

Example inputs: [testData/rcsInputs/](../testData/rcsInputs/).

```sh
mpiexec -n 1 ./build/gnu-release-mpi/bin/opensemba_rcs \
  -i testData/rcsInputs/3D_Nasa_Almond_G2_25cm_5GHz.rcs.json
```

### Top-level fields

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `runmode` | string | yes | — | Export path segment under `exports/SimulationData/` (must match how the simulation was launched, e.g. `cuda-1`, `single-core`, `mpi-8`). |
| `casename` | string | yes | — | Case folder / JSON base name under `testData/maxwellInputs/`. |
| `frequencies` | object | yes | — | Frequency linspace in **Hz** (see below). |
| `angles` | object | yes | — | Spherical angle grids in **radians** (see below). |
| `max_time` | double | no | (none) | Keep only snapshots with simulation time ≤ this value. |
| `every_n_steps` | integer | no | `1` | Keep snapshot index `i` if `i % every_n_steps == 0` while reading. Skipped payloads are not loaded. **Not** equivalent to the full time series (Nyquist / spectrum change). Must be ≥ 1. |
| `ram_gate` | double | no | (none) | Hard memory budget in **GiB**. If omitted, budget ≈ 50% of `MemAvailable`. If estimated load-all peak exceeds the budget, use streaming DFT. Use a small value to force streaming for debugging (e.g. `0.05`). Must be `> 0` when set. |

### frequencies

| Field | Type | Description |
|-------|------|-------------|
| `start` | double | First frequency (Hz). |
| `end` | double | Last frequency (Hz). |
| `steps` | integer | Number of samples (`linspace`; `1` ⇒ only `start`). |

### angles

| Field | Type | Description |
|-------|------|-------------|
| `theta` | object | `{ start, end, steps }` in radians. |
| `phi` | object | `{ start, end, steps }` in radians. |

Each of `theta` / `phi` uses the same linspace convention as `frequencies`. The post-processor evaluates the Cartesian product of all θ and φ samples.

### Minimal example

```json
{
  "runmode": "cuda-1",
  "casename": "3D_Nasa_Almond_G2_25cm_5GHz",
  "ram_gate": 100,
  "frequencies": { "start": 1e9, "end": 5e9, "steps": 3 },
  "angles": {
    "theta": { "start": 1.570796326794896, "end": 1.570796326794896, "steps": 1 },
    "phi": { "start": 0.0, "end": 6.28318530718, "steps": 361 }
  }
}
```

### Outputs

Under `exports/SimulationData/<runmode>/<casename>/RCSSurface/<probe_name>/`:

| Path | Contents |
|------|----------|
| `farfield/farfieldData_Th_*_Phi_*_dgtd.dat` | Far-field radiation potential vs frequency |
| `rcs/rcsData_Th_*_Phi_*_dgtd.dat` | RCS vs frequency |

One file pair per requested angle.
