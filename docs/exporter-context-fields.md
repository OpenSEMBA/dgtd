# Exporter context fields

The ParaView exporter (`probes.exporter`) always writes the time-domain
electric and magnetic fields `E` and `H`. With `context_fields` it also writes
additional **DOF-attached** fields that describe the model: material constants,
dispersive plasma parameters, PML stretch profiles, and mesh tags. They give
context for what the E/H fields are propagating through — for example, the
local plasma frequency `f_p` and collision frequency `nu_c` on a RAM-C II
sheath.

This page explains what each field means, how the values are attached to the
DG degrees of freedom, and how the exporter decides which fields to write.
The JSON keys are also listed in
[json-input-format.md](./json-input-format.md#exporter-context-fields).

## How the values are attached

Every context field is a scalar `mfem::ParGridFunction` on the same
discontinuous Galerkin space as `E`/`H`. For each element the exporter:

1. Reads the element attribute and the element DOF list.
2. For each local DOF, obtains the physical coordinate from the finite
   element's nodal point and the element transformation
   (`fe->GetNodes().IntPoint(j)` + `tr->Transform`).
3. Evaluates the quantity at that coordinate and stores it in the DOF.

Consequences:

- **Element-constant** quantities (`eps_r`, `mu_r`, `sigma_S_per_m`,
  `mat_tag`, `eps_inf`, `eps_s`, `tau_s`, `h_m`) are written to every DOF of
  the element.
- **Spatial profiles** (`f_p`, `f_1`, `nu_c` from `gaussian`,
  `positive_gaussian_derivative`, `gridfunction`, or `ramc_*` profiles) are
  evaluated at each DOF, so they vary DOF by DOF. This is required for plasma
  sheaths whose parameters are not element-constant. A 1D profile (for
  example `ramc_50km_*`) is stretched over the tagged region's `x` extent and
  is uniform in the transverse directions.
- PML profiles are evaluated at each DOF and for each active axis.
- The values are **nodal evaluations** at the DG points, not the
  mass-matrix-projected operator coefficients `M_K^{-1} M_K[f]`. For smooth
  profiles the two agree; for steep profiles the nodal values are the natural
  visualization and avoid projection ringing.

## JSON keys

```json
"probes": {
  "exporter": {
    "save_every": 0.1,
    "context_fields": "all",
    "context_mode": "auto"
  }
}
```

| Key | Default | Meaning |
|---|---|---|
| `context_fields` | `"all"` | `"all"`, `"none"`, or an array of field names. `"all"` selects only the fields relevant to the model (see [Applicability](#applicability)). An explicit array is honored as written. |
| `context_mode` | `"auto"` | `"auto"`, `"inline"`, `"separate"`, or `"off"`. |

### `context_mode`

- `"auto"` and `"inline"`: time-independent context fields are written in the
  **first frame** of the main `<name>.pvd` (at `t = 0`) and then removed from
  the collection, so later frames carry only `E`/`H`. This keeps the dataset
  in a **single `.pvd`**. If a future field is time-dependent, it is written
  every cycle instead.
- `"separate"`: writes the context fields to a separate
  `<name>_context` collection (one frame if time-independent, every cycle if
  time-dependent).
- `"off"`: disables context fields (same as `context_fields: "none"`).

The solver prints the decision at startup, for example:

```text
[Exporter] MaxwellView: 7 context fields, time-independent -> t=0 (inline, single .pvd)
```

## Field reference

All units are SI.

| Field | Units | Meaning | Source |
|---|---|---|---|
| `eps_r` | – | Relative permittivity (on dispersive tags this is `eps_inf`). | `Material::getPermittivity()` |
| `mu_r` | – | Relative permeability. | `Material::getPermeability()` |
| `sigma_S_per_m` | S/m | Bulk conductivity. | `Material::getConductivity() / Z0_SI` |
| `mat_tag` | – | Element attribute (mesh material tag). | `mesh.GetAttribute(el)` |
| `h_m` | m | Element size. | `mesh.GetElementSize(el)` |
| `eps_inf` | – | Instantaneous (high-frequency) relative permittivity of a Debye/Lorentz/Drude tag. | `DebyeProperties` / `LorentzProperties` |
| `fp_Hz` | Hz | Plasma frequency `f_p = omega_p / (2*pi)`. Profile-aware. | `LorentzProperties::evalOmegaP` |
| `f1_Hz` | Hz | Lorentz resonance `f_1 = omega_1 / (2*pi)`; 0 for Drude. Profile-aware. | `LorentzProperties::evalOmega1` |
| `nuc_Hz` | Hz | Collision frequency `nu_c = gamma / (2*pi)`. Profile-aware. | `LorentzProperties::evalGamma` |
| `eps_s` | – | Static relative permittivity of a Debye tag. | `DebyeProperties::eps_s` |
| `tau_s` | s | Debye relaxation time. | `DebyeProperties::tau_solver / c_SI` |
| `pml_sigma_x/y/z` | 1/m | SC-PML stretch conductivity per axis. | `PMLDirectionProfiles::sigma` |
| `pml_kappa_x/y/z` | – | SC-PML stretch factor. | `PMLDirectionProfiles::kappa` |
| `pml_alpha_x/y/z` | – | SC-PML shift. | `PMLDirectionProfiles::alpha` |
| `pml_depth_x/y/z` | m | SC-PML depth per axis. | `PMLDirectionProfiles::depth` |
| `pml_mask` | – | 1 inside a volumetric PML, 0 elsewhere. | `Model::isPMLAttribute` |

### Physical interpretation

- `eps_r`, `mu_r`, `sigma_S_per_m` are the usual bulk material constants. On a
  Debye/Lorentz/Drude tag, `eps_r` is `eps_inf`; use `eps_inf` for clarity.
- `fp_Hz` is the plasma frequency in ordinary Hz (`omega_p / 2*pi`). It is the
  frequency below which a cold plasma reflects; above it the plasma is
  transparent.
- `nuc_Hz` is the collision frequency in ordinary Hz (`gamma / 2*pi`). It sets
  the collisional absorption.
- `f1_Hz` is the Lorentz resonance in ordinary Hz. It is zero for a Drude
  (cold-plasma) pole and is omitted from `"all"` when every pole is Drude.
- `eps_s` and `tau_s` describe a single-pole Debye relaxation.
- The PML fields are the normalized Cartesian SC-PML stretch profiles used by
  the operator. `pml_sigma_*` has units of 1/m; `pml_kappa_*` and
  `pml_alpha_*` are dimensionless; `pml_depth_*` is in metres. Inactive axes
  are not written when `context_fields: "all"` is used.

## Smooth plasma visualization

When `probes.exporter.smooth_plasma_fields` is `true`, the DG rates are
complemented by two visualization-only collections.  The solver
still reads the original `omega_p.gf` / `nu_c.gf` values at the DG quadrature
points; neither collection feeds back into the operator.

- `<name>_plasma` is a continuous H1 order-2 reconstruction of `fp_Hz` and
  `nuc_Hz`.  The DG rate is sampled at the H1 nodes, so the sub-element
  variation is retained; shared nodes are volume-averaged across the incident
  plasma elements.  The collection is written on a **plasma-only submesh**
  (the Lorentz tags), so it cannot extend into vacuum.  The continuous field
  is floored at the minimum element-centre rate over the plasma volume (the
  free-stream value), so DG polynomial undershoots cannot create artificial
  zero or negative spots in ParaView.  This is a visualization floor, not a
  change to the solver rates.
- `<name>_plasma_cells` is the piecewise-constant view on the same submesh.
  `fp_Hz_cell` / `nuc_Hz_cell` are the DG rate at the element centroid (a
  local representative value), and `fp_Hz_mean` / `nuc_Hz_mean` are the exact
  L2 element means used by the DG mass integration.  The element-centre
  sample is the better default for a coarse cell that spans a strong gradient;
  the mean is available for solver-side comparisons.

A continuous view is intentionally preferred for the sheath: an element mean
is a volume average, and a single mesh element can span a large part of the
sheath (the RAM-C II nose stagnation element is a 24 mm axial sliver).  The
mean then assigns a high value to the low-density nose vertex and makes the
cell view look blocky even though the underlying DG field carries the local
variation.

## Applicability

`context_fields: "all"` (the default) writes only fields that exist in the
model:

- `eps_r`, `mu_r`, `mat_tag`, `h_m` are always written.
- `sigma_S_per_m` is written only when at least one material has a nonzero
  bulk conductivity.
- `eps_inf`, `fp_Hz`, `nuc_Hz` are written when the model has a Lorentz or
  Drude pole.
- `f1_Hz` is written only when at least one Lorentz pole has a nonzero
  `omega_1` (or an `omega_1` spatial profile). Pure Drude models omit it.
- `eps_s`, `tau_s` are written when the model has a Debye pole.
- `pml_sigma_*`, `pml_kappa_*`, `pml_alpha_*`, `pml_depth_*`, `pml_mask` are
  written only when the model has a volumetric PML, and only for the active
  stretch axes.

An explicit array such as
`"context_fields": ["eps_r", "mu_r", "fp_Hz", "nuc_Hz"]` is honored as
written; fields that are not applicable to the model keep their default
values (0, or 1 for `mu_r`/`eps_r`/`pml_kappa_*`).

## Examples

### Lorentz plasma with a spatial profile

```json
"probes": {
  "exporter": {
    "save_every": 0.05,
    "context_fields": "all",
    "context_mode": "auto"
  }
}
```

With a Lorentz pole present, `"all"` writes `eps_r`, `mu_r`, `mat_tag`,
`h_m`, `eps_inf`, `fp_Hz`, `f1_Hz`, `nuc_Hz`. If `omega_p`/`gamma` are
`gridfunction` profiles, `fp_Hz`/`nuc_Hz` vary DOF by DOF.

### Drude (cold plasma)

A `drude` material has `omega_1 = 0`, so `f1_Hz` is not written. The fields
are `eps_r`, `mu_r`, `mat_tag`, `h_m`, `eps_inf`, `fp_Hz`, `nuc_Hz`.

### Debye

A Debye material writes `eps_r`, `mu_r`, `mat_tag`, `h_m`, `eps_inf`,
`eps_s`, `tau_s`.

### PML

A volumetric PML adds the `pml_*` fields for its active axes and
`pml_mask`. For a 1D PML only `pml_sigma_x`, `pml_kappa_x`, `pml_alpha_x`,
`pml_depth_x`, and `pml_mask` are written.

### Explicit subset

```json
"context_fields": ["eps_r", "mu_r", "fp_Hz", "nuc_Hz"]
```

### Disable

```json
"context_mode": "off"
```

or

```json
"context_fields": "none"
```

## Units and the `gamma`/`nu_c` gridfunction

The solver stores dispersive rates in normalized units (SI rad/s divided by
`c_SI`). The exporter converts them back to ordinary SI:

- `fp_Hz = omega_p_solver * c_SI / (2*pi)`
- `f1_Hz = omega_1_solver * c_SI / (2*pi)`
- `nuc_Hz = gamma_solver * c_SI / (2*pi)`
- `tau_s = tau_solver / c_SI`
- `sigma_S_per_m = sigma_solver / Z0_SI`

For `gridfunction` rate profiles, the file values are treated as **rad/s** by
default. If the file stores ordinary frequencies in Hz, add `"unit": "Hz"` to
the profile object; the solver multiplies by `2*pi` before the `1/c_SI`
conversion. This is relevant for `nu_c` archives written in 1/s: use
`"unit": "Hz"` so `nuc_Hz` reproduces the file's collision frequency.

## Practical notes

- The context fields are static in the current implementation. The single
  `.pvd` therefore has them only in the first frame; step through the time
  series for `E`/`H` and inspect the first frame for the material context.
- ParaView can contour any scalar field directly. `mat_tag` is a useful
  selector for clipping to a material region.
- The TFSF interface is a total-field/scattered-field discontinuity, not a
  context-field artifact. The total-field region is smooth when the mesh
  resolves the incident pulse.
- For a `"separate"` collection, the files are written under
  `exports/ParaView/<run-mode>/<name>_context/`.
