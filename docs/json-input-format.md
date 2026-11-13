# JSON input format

Simulation cases are defined by JSON files parsed at runtime. Each case must follow the repository naming layout: the case folder name, the `.json` filename, the `.msh` filename, and `model.filename` must all share the same `<case_name>` base, and `model.filename` must not include a path.

Example: [testData/testInputs/1D_PEC/1D_PEC.json](../testData/testInputs/1D_PEC/1D_PEC.json)

Legend: **[REQUIRED]** = must be present; **[OPTIONAL]** = has a default value (shown).

## solver_options

Object. User can customise solver settings. If undefined, all defaults apply.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `evolution_operator` | string | `"global"` | `"global"` (default; assembled CSR, wires, implicit ODE). `"hesthaven"` is the element-local operator for conductivity, volumetric PML, Debye, Lorentz, SGBC, interior boundaries, and centered SMA (2D/3D). `"maxwell"` is deprecated. |
| `upwind_alpha` | double | `1.0` | Upwind flux blending: `0.0` = centered, `1.0` = upwind. |
| `final_time` | double | `2.0` | Simulation duration in natural units (1 m/c). |
| `time_step` | double | `0.0` | Fixed time step in natural units. Required for 2D/3D. In 1D, `0.0` triggers automatic CFL-based step. |
| `cfl` | double | `1.0` | CFL for automatic 1D time step. Range (0, 1]. Ignored if `time_step` is set. Not used in 2D/3D auto step. |
| `sgbc_cfl` | double | `0.5` | Crossing-time CFL for SGBC sub-step recommendation: `δt_rec = crossing_time × sgbc_cfl × opacity_relax` (then SI→natural), where `opacity_relax = min((d/δ)^2, 50)` for thick conductive layers ($d/\delta > 3$ evaluated at the source $1/e$ power frequency). Must be `> 0`. Omit for historical default (`0.5`). Does not change SGBC layer geometry; only `nsteps = ceil(Δt / δt_rec)`. |
| `order` | integer | `2` | Polynomial order of the FE basis. |
| `spectral` | boolean | `false` | Spectral evolution operator (full matrix eigenvalue step). High cost; limited feature support. |
| `export_operator` | boolean | `false` | Write the assembled evolution operator to `exports/Operators/<mesh>/` (`{mesh}_global.csr`, and `{mesh}_tfsf.csr`, `{mesh}_tfsf_mapping.csr`, `{mesh}_farfield.csr` when those pieces exist). With more than one MPI rank the partitioned matrices are gathered into the same serial-mesh DOF order a one-rank run writes. The directory name does not include the run mode. |
| `basis_type` | integer | `1` | MFEM basis: `0` GaussLegendre, `1` GaussLobatto, `2` Bernstein, `3` OpenUniform, `4` CloseUniform, `5` OpenHalfUniform. |
| `ode_type` | integer | `0` | Time integrator: `0` RK4, `1` BackwardEuler, `2` Trapezoidal, `3` ImplicitMidpoint, `4` SDIRK33, `5` SDIRK23, `6` SDIRK34. |
| `checkpoint_percent` | double | `20` | `0` disables checkpoints. Otherwise, after the completed step that crosses each multiple of this percent of `final_time`, write the evolved state. The default `20` saves near 20%, 40%, 60%, 80%, and 100%. |

### Checkpoints

`checkpoint_percent` writes under `exports/SimulationData/<run-mode>/<case>/Checkpoints/`. Each save stores the packed field state (E, H, and any PML, Debye, or Lorentz auxiliaries), SGBC face history, the time and time step, probe cursors, the element partition, and the solver run time up to that save. Operators are rebuilt on the next launch. A finished run deletes that directory. Probe files, ParaView output, and `SimulationStats` stay.

`Simulation Run Time` in `statistics_rank*.dat` is the time spent inside the time loop and covered by a checkpoint, plus the time after the last resume. Minutes after the last save and before a crash are not included. A run that reached its 40% save in 4 minutes, died at 50%, and then took 6 minutes to finish after `--continue` reports 10 minutes.

Resume with the same JSON, the same mesh, and the same MPI rank count and device, from the same working directory:

```sh
mpirun -np N ./opensemba_dgtd -i case.json --continue
```

`solver_options.checkpoint_percent` may differ on resume. Any other JSON change, or a different mesh, is rejected. A checkpoint written with another `-np` or device is rejected as well, including when `--continue` is omitted, so that a plain launch does not start at t = 0. `--restart` deletes every `Checkpoints` directory for this case and starts from t = 0:

```sh
mpirun -np N ./opensemba_dgtd -i case.json --restart
```

`SIGINT` or `SIGTERM` delivered to the `opensemba_dgtd` processes finishes the current step, writes one checkpoint, and exits. `mpirun` handles Ctrl-C itself and can kill those processes before that write; the last complete checkpoint is what a later `--continue` loads. A failed save does not stop the run. The previous complete checkpoint is kept, and the next attempt waits at least a minute.

### evolution_operator: hesthaven

`"hesthaven"` is element-local dense operators (matrix-free on straight meshes). It runs piecewise-constant materials, bulk conductivity, Cartesian volumetric PML, Debye, and Lorentz/Drude, including spatial rate profiles, as per-element updates. Interior-boundary fluxes and centered SMA (`upwind_alpha: 0`) use the assembled face operator on the elements that touch those faces. SGBC reuses the same implicit coating sub-step as Global and still requires RK4. It does **not** run thin wires, delta-gap, or coaxial ports. Other shared explicit physics can use an implicit `ode_type`. It is 2D/3D only.

| Capability | hesthaven | global |
|------------|-----------|--------|
| MPI (`mpirun`) | Yes — shared-face ghost exchange and neighbor connectivity in `Mult()` | Yes |
| CUDA | Yes — CUDA builds of `opensemba_dgtd` default to `--device cuda`; pass `--device cpu` (or `omp`) to stay on the host | Yes |
| Conductivity / PML / Debye / Lorentz | Yes — element-local, 2D/3D | Yes |
| SGBC | Yes — same implicit sub-step, parent RK4 | Yes |
| Thin wires | No | Yes |
| Implicit `ode_type` | Yes, except SGBC | Yes, except SGBC |
| Centered SMA (`upwind_alpha: 0`) | Yes | Yes |

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
| `debye` | object | — | Single-pole electric Debye on an untyped material. See below and [dispersive-materials.md](./dispersive-materials.md). |
| `lorentz` | object | — | Single-pole electric Lorentz on an untyped material. See below and [dispersive-materials.md](./dispersive-materials.md). |
| `drude` | object | — | Cold-plasma Drude pole ($\omega_1=0$) on an untyped material. See below and [dispersive-materials.md](./dispersive-materials.md). |

#### debye

Optional object. Legal only when the material entry has no `type`. The electric mass uses `eps_inf`. `tau` is seconds; the solver stores $\tau_{\mathrm{SI}}\,c_{\mathrm{SI}}$ because time is in light-meters.

| Field | Type | Description |
|-------|------|-------------|
| `eps_inf` | double | Instantaneous relative permittivity. Must be $\ge 1$. |
| `eps_s` | double | Static relative permittivity. Must be greater than `eps_inf`. |
| `tau` | double | Relaxation time in seconds. Must be $> 0$. |

`relative_permittivity` is rejected together with `debye`. `relative_permeability` and `bulk_conductivity` keep their usual defaults. `bulk_conductivity` stays an independent Ohm term beside the pole.

`debye` is rejected on `type: "vacuum"`, `type: "PML"`, and on an SGBC layer. Initialization also aborts if a Debye tag is a PML tag. A tag cannot be listed as both Debye and another material. `evolution_operator` must be `"global"` or `"hesthaven"` (`"hesthaven"` is 2D/3D). `spectral: true` is rejected. Implicit ODE includes the Debye polarization.

#### lorentz

Optional object. Legal only when the material entry has no `type`. The electric mass uses `eps_inf`. Give the three rates in **one** style, never mixed. The solver stores $\omega/c_{\mathrm{SI}}$ because time is in light-meters. $f_1 = 0$ (or $\omega_1 = 0$) is a cold plasma on this same pole.

The relative permittivity is

$$
\varepsilon_r(\omega)=\varepsilon_\infty+\frac{\omega_p^2}{\omega_1^2-\omega^2+j\,\gamma\omega}.
$$

**Hz style (preferred for new cases):** `eps_inf`, `f_p`, `f_1`, `gamma`. All three rates are ordinary frequencies in Hz. Conversion is $\omega=2\pi f$, including damping: JSON `gamma: 2.5e7` is $25\,\mathrm{MHz}$, stored as $\gamma=2\pi\times 2.5\times 10^7\,\mathrm{rad/s}$.

**rad/s style:** `eps_inf`, `omega_p`, `omega_1`, `gamma`. All three rates are angular frequencies in rad/s. Do not mix with `f_p` / `f_1`.

| Field | Type | Description |
|-------|------|-------------|
| `eps_inf` | double | Instantaneous relative permittivity. Must be $\ge 1$. |
| `f_p` | double or object | Plasma frequency in Hz (or spatial profile object). Must be $> 0$. |
| `f_1` | double or object | Resonance frequency in Hz (or spatial profile object). Must be $\ge 0$. |
| `omega_p` | double or object | Plasma frequency in rad/s (or spatial profile object). Must be $> 0$. |
| `omega_1` | double or object | Resonance frequency in rad/s (or spatial profile object). Must be $\ge 0$. |
| `gamma` / `nu_c` | double or object | Collision frequency in $\ddot{P}+\gamma\dot{P}+\omega_1^2 P=\omega_p^2 E$. Hz with `f_*`, or rad/s ($1/\mathrm{s}$) with `omega_*` (or spatial profile object). Must be $\ge 0$. |

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

`relative_permittivity` is rejected together with `lorentz`. `bulk_conductivity` stays an independent Ohm term. A tag cannot carry both `debye` and `lorentz` (or `drude`). `lorentz` is rejected on vacuum, PML, and SGBC. Initialization also aborts if a Lorentz tag is a PML tag. `evolution_operator` must be `"global"` or `"hesthaven"` (`"hesthaven"` is 2D/3D). `spectral: true` is rejected. Implicit ODE includes the Lorentz auxiliaries.

#### drude

Optional object. Shorthand for a cold-plasma Lorentz pole with $\omega_1 = 0$:

$$
\varepsilon_r(\omega)=\varepsilon_\infty-\frac{\omega_p^2}{\omega^2-j\,\nu_c\omega}.
$$

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `eps_inf` | double | `1.0` | Instantaneous relative permittivity ($\ge 1$). |
| `omega_p` / `f_p` | double or object | — | Plasma frequency in rad/s (`omega_p`) or Hz (`f_p`), or spatial profile object. |
| `nu_c` / `gamma` | double or object | — | Collision frequency in $1/\mathrm{s}$ (with `omega_p`) or Hz (with `f_p`), or spatial profile object. |

#### Spatial rate profiles (`lorentz` and `drude`)

Any rate in `lorentz` or `drude` (`omega_p`, `f_p`, `omega_1`, `f_1`, `gamma`, `nu_c`) can be either a constant number or a spatial profile object `{ "type": "...", ... }`. When a profile object is used on an attribute, the solver evaluates the profile at each element's quadrature points via `mfem::FunctionCoefficient` inside `mfem::MassIntegrator` and pre-multiplies by the element inverse mass $M_K^{-1}$. `GlobalEvolution` stores those blocks in `lorentzOperator_`. `HesthavenEvolution` stores the same blocks per element.

##### 1D automatic region stretching

In 1D, the solver scans the mesh at initialization for all elements carrying the material's `tags` and records the physical bounding interval $[x_{\min}^{\mathrm{mesh}}, x_{\max}^{\mathrm{mesh}}]$. Every quadrature point $x \in [x_{\min}^{\mathrm{mesh}}, x_{\max}^{\mathrm{mesh}}]$ is mapped to the normalized unit coordinate:

$$
\xi(x) = \frac{x - x_{\min}^{\mathrm{mesh}}}{x_{\max}^{\mathrm{mesh}} - x_{\min}^{\mathrm{mesh}}} \in [0, 1].
$$

Each profile defines its canonical shape on $\xi \in [0, 1]$, so the complete profile automatically stretches or compresses to fit the tagged mesh region whether that region is $10\,\mathrm{mm}$, $12\,\mathrm{mm}$, $24\,\mathrm{mm}$, or $2\,\mathrm{m}$ wide.

##### Supported profile `type` values

| `type` | Required keys | Optional keys (defaults) | Description |
|--------|---------------|--------------------------|-------------|
| `"constant"` | `value` | — | Uniform rate across the region. |
| `"gaussian"` | `peak` (or `amplitude`) | `cutoff_factor` (`1e-2`), `offset` (`0.0`) | Symmetric bell curve stretched between its `cutoff_factor` tails. |
| `"positive_gaussian_derivative"` | `peak` (or `amplitude`) | `cutoff_factor` (`1e-2`), `offset` (`0.0`) | Rectified two-lobe Gaussian derivative stretched between its outer `cutoff_factor` tails. |
| `"gridfunction"` | `mesh`, `file` | `axisymmetric` (`true`), `unit` (`"rad/s"`), `fallback` (`"nearest"`), `default` (`0.0`) | Scalar H1 GridFunction from a plasma_kit archive (`mesh.mesh` + `.gf`), evaluated pointwise at the DG quadrature points. |
| `"ramc_61km_omega_p"` | — | `direction` (`1.0`) | RAM-C II 61 km stagnation-line $\omega_p(x)$ $[\mathrm{rad/s}]$ mapped from $[-0.012, 0.0]\,\mathrm{m}$ onto $\xi \in [0, 1]$. |
| `"ramc_61km_nu_c_gupta"` | — | `direction` (`1.0`) | RAM-C II 61 km stagnation-line $\nu_c(x)$ $[1/\mathrm{s}]$ (Gupta closure) mapped onto $\xi \in [0, 1]$. |
| `"ramc_61km_nu_c_liu"` | — | `direction` (`1.0`) | RAM-C II 61 km stagnation-line $\nu_c(x)$ $[1/\mathrm{s}]$ (Liu closure) mapped onto $\xi \in [0, 1]$. |
| `"ramc_50km_omega_p"` | — | `direction` (`1.0`) | RAM-C II 50 km stagnation-line $\omega_p(x)$ $[\mathrm{rad/s}]$ mapped from $[-0.024, 0.0]\,\mathrm{m}$ onto $\xi \in [0, 1]$. |
| `"ramc_50km_nu_c_gupta"` | — | `direction` (`1.0`) | RAM-C II 50 km stagnation-line $\nu_c(x)$ $[1/\mathrm{s}]$ (Gupta closure) mapped onto $\xi \in [0, 1]$. |
| `"ramc_50km_nu_c_liu"` | — | `direction` (`1.0`) | RAM-C II 50 km stagnation-line $\nu_c(x)$ $[1/\mathrm{s}]$ (Liu closure) mapped onto $\xi \in [0, 1]$. |

##### `"gridfunction"` profiles

A `"gridfunction"` profile reads an MFEM mesh and a scalar H1 `GridFunction`
written by plasma_kit (`omega_p.gf`, `nu_c.gf`). The file must store SI
`rad/s` (or `1/s`); `unit: "Hz"` multiplies by `2*pi`. Only the `1/c0`
conversion to solver units is applied. Relative `mesh`/`file` paths are
resolved next to the case JSON first and then from the current directory.

```json
"lorentz": {
  "eps_inf": 1.0,
  "f_p": {
    "type": "gridfunction",
    "mesh": "ramc2_subghz_2d_gupta/mesh.mesh",
    "file": "ramc2_subghz_2d_gupta/omega_p.gf",
    "axisymmetric": true
  },
  "f_1": 0.0,
  "gamma": {
    "type": "gridfunction",
    "mesh": "ramc2_subghz_2d_gupta/mesh.mesh",
    "file": "ramc2_subghz_2d_gupta/nu_c.gf",
    "axisymmetric": true
  }
}
```

The rate profiles are evaluated at every element quadrature point and passed
through the same `M_K^{-1} M_K[f]` projection as the built-in profiles, so the
H1 data is effectively projected onto the DG (L2) space of the simulation,
element by element, only on the Lorentz material's volume tags. For a 2D
source mesh, the 3D point `(x, y, z)` is reduced to `(x, r = hypot(y, z))`
when `axisymmetric` is true (the default). Points outside the source region
return the nearest source element value (`fallback: "nearest"`, default) or
`default` (0) when `fallback: "default"`.

##### How `"gaussian"` and `"positive_gaussian_derivative"` work

Because plasma, resonance, and collision frequencies ($\omega_p, \omega_1, \nu_c$) must be non-negative and self-contained on arbitrary mesh intervals, both profiles use `cutoff_factor` $\epsilon \in (0, 1)$ (default `1e-2`) to determine their canonical span $[-z_{\mathrm{cut}}, +z_{\mathrm{cut}}]$ and map it linearly onto $\xi \in [0, 1]$ via $z(\xi) = (2\xi - 1)\,z_{\mathrm{cut}}$:

1. **`"gaussian"`** evaluates the canonical Gaussian $G(z) = \exp(-z^2 / 2)$:

$$
z_{\mathrm{cut}} = \sqrt{-2\ln(\epsilon)} \quad (\approx 3.035\text{ for }\epsilon = 10^{-2}),
$$

$$
f(\xi) = \mathrm{offset} + \mathrm{peak}\cdot\exp\!\left(-\frac{z(\xi)^2}{2}\right).
$$

At the region boundaries $\xi = 0$ and $\xi = 1$, $f$ equals $\mathrm{offset} + \epsilon\cdot\mathrm{peak}$; at the region midpoint $\xi = 0.5$ ($z = 0$), $f$ reaches $\mathrm{offset} + \mathrm{peak}$.

2. **`"positive_gaussian_derivative"`** takes the magnitude of the normalized Gaussian derivative so that both lobes are positive and have unit peak height at $z = \pm 1$:

$$
D(z) = |z|\,\exp\!\left(\frac{1 - z^2}{2}\right), \qquad f(\xi) = \mathrm{offset} + \mathrm{peak}\cdot D\big(z(\xi)\big).
$$

The outer cutoff $z_{\mathrm{cut}} > 1$ is solved at initialization from:

$$
z_{\mathrm{cut}}\,\exp\!\left(\frac{1 - z_{\mathrm{cut}}^2}{2}\right) = \epsilon \quad (\approx 3.648\text{ for }\epsilon = 10^{-2}).
$$

Across the normalized region $\xi \in [0, 1]$, `"positive_gaussian_derivative"` therefore produces:
- $\xi = 0$ ($x = x_{\min}^{\mathrm{mesh}}$, $z = -z_{\mathrm{cut}}$): outer tail $\mathrm{offset} + \epsilon\cdot\mathrm{peak}$
- $\xi = \tfrac{1}{2}(1 - 1/z_{\mathrm{cut}}) \approx 0.363$ ($z = -1$): first lobe peak $\mathrm{offset} + \mathrm{peak}$
- $\xi = 0.5$ ($x = \text{midpoint}$, $z = 0$): central trough $\mathrm{offset}$
- $\xi = \tfrac{1}{2}(1 + 1/z_{\mathrm{cut}}) \approx 0.637$ ($z = +1$): second lobe peak $\mathrm{offset} + \mathrm{peak}$
- $\xi = 1$ ($x = x_{\max}^{\mathrm{mesh}}$, $z = +z_{\mathrm{cut}}$): outer tail $\mathrm{offset} + \epsilon\cdot\mathrm{peak}$

Examples:

```json
"drude": {
  "eps_inf": 1.0,
  "omega_p": {
    "type": "positive_gaussian_derivative",
    "peak": 4.0e9,
    "cutoff_factor": 1.0e-2
  },
  "nu_c": {
    "type": "gaussian",
    "peak": 1.0e9,
    "cutoff_factor": 1.0e-2
  }
}
```

RAM-C II 61 km example ([1D_Drude_RAMC_61km](../testData/testInputs/1D_Drude_RAMC_61km/1D_Drude_RAMC_61km.json)):

```json
"drude": {
  "eps_inf": 1.0,
  "omega_p": {
    "type": "ramc_61km_omega_p"
  },
  "nu_c": {
    "type": "ramc_61km_nu_c_gupta"
  }
}
```

### wires [OPTIONAL]

Array. Each entry is one Holland thin wire. The curve is a linear Gmsh line (element type 1) saved in the volume `.msh`. Its elementary tag is the JSON tag. The solver reads those lines from the case file, then builds the volume mesh from a temporary copy with the lines removed. The `.msh` you launch is left unchanged.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tags` | int[] | — | Elementary tags of the line elements that form this curve. |
| `radius` | double | — | Wire radius in metres. Must be positive and smaller than $\Delta$ (`transverseStep`, or the host mean edge when that is omitted). |
| `resistancePerMeter` | double | `0.0` | Extra series resistance, ohm/m. |
| `inductancePerMeter` | double | `0.0` | Extra series inductance, H/m, added to Edelvik's in-cell inductance. |
| `transverseStep` | double | host mean edge | Transverse length $\Delta$ in metres. The tube radius is $1.7\Delta$ and the in-cell inductance is $(1/(2\pi))\log((1.7\Delta+a)/(2a))$. Omit it to use the host tetrahedron's mean edge. |

A `wire_voltage` source, when present, adds $V(t)/\ell$ to the tangential electric drive of the one piece that contains `position`. $\ell$ is that piece's length. $V$ is the Gaussian below, in the same units as $\int E\,\mathrm{d}\ell$.

The wire stays in vacuum or in a PML region whose `kappa_max` is `1`. With that value the stretch factor on `E` is 1, so the wire reads the physical tangential field and returns its current on the electric rows after the PML update. `kappa_max` above 1 is rejected. A dispersive material, or an element that touches TFSF, SGBC, or a delta-gap, is rejected. A free end sets the current to zero. An end on PEC sets the charge per length to zero. An end on PMC, SMA, or SGBC is rejected. The piece chain and the line current are replicated on every MPI rank. Each rank adds the tangential electric field of the tetrahedra it owns, and only that rank adds the return current on those elements. A CUDA run is one MPI rank. The line update then stays on the device and does not copy the volume state. `ode_type` 0 (RK4) is required. `evolution_operator` must be `"global"`. A mid-wire current time series is `probes.wire`, written by rank 0. That current is the light-metre unknown, not amperes. Dividing it by $Z_0$ is a label only.

When `probes.exporter` is set, each wire is written inside that volume case folder. Next to `3D_Wire_Box.msh/3D_Wire_Box.msh.pvd` the run adds `3D_Wire_Box.msh/3D_Wire_Box_wire00.msh/3D_Wire_Box_wire00.msh.pvd`, then `wire01`, and so on. Each of those is the 1D mesh of that wire, with nodal fields `current` and `charge`. Open it in the same session as the volume file.

### boundaries [REQUIRED]

Array. At least one entry per boundary or interface condition.

| Field | Type | Description |
|-------|------|-------------|
| `tags` | int[] | Boundary attribute IDs. |
| `type` | string | `"PEC"`, `"PMC"`, `"SMA"`, or `"SGBC"`. |

#### type: SGBC

Surface general boundary condition (multi-layer absorber / material slab sub-solver). Full guide (automatic layer meshing, Péclet/skin-depth rules, sub-stepping, and ParaView export): [sgbc.md](./sgbc.md).

| Field | Description |
|-------|-------------|
| `exporter_probe` | Legacy optional flag (accepted for backward compatibility). When top-level `probes.exporter` exists, internal SGBC fields are automatically exported per tag as `<stem>_sgbc_tag<tag>.msh` inside the exporter directory. |
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

Each SGBC sub-solver caches up to two dense LU factorizations of $J = I - \gamma\,\delta t\,A_s$, covering both the RK4 half-stage sub-step ($\delta t_{1/2} = 0.5\Delta t / N_{1/2}$) and the full-step sub-step ($\delta t_1 = \Delta t / N_1$) without re-factorizing during the time loop. Multiple SGBC boundary blocks with distinct `tags` and layer stacks may be defined in the same model and distributed across MPI ranks.

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

### wire

Array. Current on a Holland thin wire, at a point that lies on the curve. This is the mid-wire current recorded by the FDTD `holland1981` case. The plane wave that drives that case is `sources` type `planewave`. The wire must sit in the total-field region and must not touch the total-field/scattered-field faces.

| Field | Description |
|-------|-------------|
| `position` | A point on the wire. |
| `field` | `"current"` or `"charge"`. Omitted means both are written. `"voltage"` is rejected: a plain Holland wire has no reference conductor. |
| `steps` | Record every N time steps (exclusive with `saves`) |
| `saves` | Total samples over the run |

The file is `exports/SimulationData/<run-mode>/<case>/WireProbes/WireProbe<id>.dat`. Columns are time in seconds, current, and charge per length.

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

Surface E/H snapshots for offline RCS (`surface_data.bin` with geometry header + time blocks). Full workflow and Huygens surface placement guide: [rcs-and-farfield.md](./rcs-and-farfield.md).

| Field | Description |
|-------|-------------|
| `tags` | Surface integration tags |
| `name` | Default `"RCSSurfaceProbe"` |
| `steps` / `saves` | As above. Prefer matching temporal density across cases (e.g. with `time_step` 0.0001 use `steps: 5` to match a 0.0005 / `steps: 1` export). |

Offline frequency/angle sweeps use a **separate** JSON for `opensemba_rcs` — see [Offline RCS JSON](#offline-rcs-json-opensemba_rcs). Do not confuse offline `every_n_steps` with case-JSON `probes.rcssurface.steps` (export cadence during `opensemba_dgtd`).

### mor_state

Full DG state snapshots in serial-mesh DOF order, the same layout as `{mesh}_global.csr`, for offline `y = A x`. With more than one MPI rank, each rank's owned degrees of freedom are gathered into that order before the file is written. Files go under `exports/SimulationData/<run-mode>/<case>/MORStateProbes/<name>/`, where `<run-mode>` is `single-core` or `mpi-N` (and `cuda-N` on a GPU run).

When the case has a planewave TFSF source, each save also writes `u_k` on the serial TFSF submesh. An MPI run builds that incident vector from the full mesh, so it matches a one-rank `u_k`.

| Field | Description |
|-------|-------------|
| `record_time_start` | Recording start time |
| `record_time_final` | Recording end time |
| `saves` | Number of uniform snapshots |
| `name` | Default `"MORState"` |

Snapshot file format (`x_0`, `x_1`, …): line 1 = time, line 2 = size, then DOF values in `[Ex | Ey | Ez | Hx | Hy | Hz]` order. `u_k` uses the same header. Its size is 6 times the number of DOFs on the TFSF submesh.

## sources [REQUIRED]

Array. At least one source; all entries superimpose. Full guide on TFSF mesh partitioning, probe placement, and source waveforms: [tfsf-and-sources.md](./tfsf-and-sources.md).

| `type` | Description |
|--------|-------------|
| `"initial"` | Volumetric initial condition |
| `"planewave"` | TFSF plane wave |
| `"dipole"` | TFSF dipole |
| `"delta_gap"` | Impressed tangential E on an interior face (2D curve or 3D surface) |
| `"coaxial_port"` | TFSF coaxial TEM mode on an interior annular face |
| `"wire_voltage"` | Series voltage on one Holland thin-wire piece |

### Gaussian width (`spread` / `f_1e`)

Time-domain Gaussians are $\exp(-(t-t_0)^2/(2\sigma^2))$ with $\sigma$ in light-metres. Set the width with **either**:

- `spread`: $\sigma$ directly, or
- `f_1e`: frequency in Hz where incident **power** is $1/e$ of DC. Then $\sigma = c/(2\pi f_{1e})$.

Do not mix units: $f_{1e}=3\times 10^8$ is 300 MHz, which is $\sigma\approx 0.16$, not `spread: 0.15` written as a frequency. If both keys are present, `f_1e` is used and a warning is printed. At solver start that 1/e frequency is compared to mean element size and FE order; if the pulse is likely under-resolved a warning is printed and the run continues.

On `planewave`, `dipole`, and `initial`, the keys live under `magnitude`. On `delta_gap`, `coaxial_port`, and `wire_voltage`, they are on the source object.

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
| `magnitude.mean` | Optional 3-vector pulse center. Omit it to place the pulse upstream of the TFSF surface. |
| `magnitude.frequency` | Optional carrier (Hz) for modulated Gaussian |
| `magnitude.file` | Sampled excitation. Two columns, time in seconds and amplitude, as in an FDTD `magnitudeFile`. Linear interpolation. The value is zero before the first sample and after the last. Time zero is the waveform on the plane $\mathbf{x}\cdot\hat{\mathbf{d}}=0$. Do not combine `file` with `spread`, `f_1e`, `frequency`, or `mean`. A relative path is read beside the JSON file. |

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

### type: wire_voltage

One entry. The mesh must already contain `model.wires`, and `evolution_operator` must be `"global"` (or omitted). `position` is a point on the wire. The piece that contains it receives the voltage. A point shared by two pieces receives half of the line integral on each piece.

`magnitude` (default `1.0`) is the peak of $V(t)$, the line integral of the extra tangential field. The field added on that piece is $V(t)/\ell$, with $\ell$ the piece length. It is held constant over the wire substeps of a stage, the same way the volume electric field is held. The pulse is $e^{-(t-t_0)^2/(2\sigma^2)}$, centered five widths after $t=0$. Width is `spread` (light-metres) or `f_1e` (Hz at $1/e$ incident power). If both are set, `f_1e` wins and a warning is printed. `signal` `"gaussian"` is that voltage. `"gaussian_derivative"` is its time derivative.

| Field | Description |
|-------|-------------|
| `position` | Required 3-vector on the wire |
| `magnitude` | Peak voltage. Default `1.0`. Must be `> 0` |
| `spread` | Gaussian width $\sigma$ in light-metres. Must be `> 0` |
| `f_1e` | Alternative to `spread`: Hz at 1/e incident power |
| `signal` | `"gaussian"` (default) or `"gaussian_derivative"` |

---

## Offline RCS JSON (`opensemba_rcs`)

Separate input file for the offline RCS / far-field post-processor (`opensemba_rcs`). It does **not** replace the Maxwell case JSON; it points at an existing export and at `testData/testInputs/<casename>/<casename>.json`, or `testData/maxwellInputs/<casename>/<casename>.json` for a local exploration case. See also [rcs-and-farfield.md](./rcs-and-farfield.md).

**Prerequisites**

1. Case JSON includes `probes.rcssurface` (see [rcssurface](#rcssurface)).
2. `opensemba_dgtd` has been run so that  
   `exports/SimulationData/<runmode>/<casename>/RCSSurface/<probe_name>/rank*/surface_data.bin`  
   (and `mesh`) exist.
3. Launch from the repository root (paths are relative: `./exports/...`, `./testData/testInputs/...`).

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
| `casename` | string | yes | — | Case folder / JSON base name under `testData/testInputs/`, or under local `testData/maxwellInputs/`. |
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
