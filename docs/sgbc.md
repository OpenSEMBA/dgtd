# Surface general boundary condition (SGBC)

Multi-layer thin-sheet and radar-absorbing material (RAM) coating boundary condition on `GlobalEvolution`.

SGBC replaces a thin conductive or multi-layer dielectric/magnetic slab on a boundary or interior surface with an embedded 1D normal-incidence DG sub-solver per face. This avoids meshing sub-millimetre layers and micron-scale skin depths inside a 2D or 3D volume mesh, and steps the stiff thin-layer equations implicitly so they do not restrict the parent explicit RK4 time step.

Use `"evolution_operator": "global"` (or omit `evolution_operator`) and `"ode_type": 0` (RK4) on the parent solver.

---

## How SGBC couples to the parent mesh

Each face tagged with `"type": "SGBC"` in `model.boundaries` instantiates a local 1D sub-mesh along the sheet normal:

```text
[ 2 left vacuum ghost cells ] | [ Layer 0 ] ... [ Layer N-1 ] | [ 2 right vacuum ghost cells ]
                              ^                               ^
                     Interior tag 3 (left)           Interior tag 4 (right)
```

At every parent RK4 stage:

1. **Parent $\to$ SGBC ghost cells:** The parent solver extracts the tangential electric and magnetic fields $(\mathbf{E}_t, \mathbf{H}_t)$ on the two sides of the SGBC face (or one side for an exterior boundary face), rotates them into a local orthonormal frame $(\hat{\mathbf{e}}_1, \hat{\mathbf{e}}_2, \hat{\mathbf{n}})$, and loads them into the left and right vacuum ghost segments.
2. **Implicit 1D sub-step (`SDIRK33`):** The 1D sub-solver advances the internal layer fields across the parent stage interval using an L-stable 3-stage SDIRK33 integrator (`GlobalEvolution::ImplicitSolve`). Because the parent uses a 4-stage explicit RK4 integrator, the SGBC solver checkpoints its internal state at the start of the time step and rolls back to that checkpoint at the beginning of each RK4 stage before advancing to the stage time.
3. **SGBC $\to$ Parent face flux:** The updated tangential traces at the left and right interfaces of the 1D stack are rotated back to global coordinates and injected into the parent RHS via the SGBC source-face operators.

### Passthrough sheets vs. backed coatings (`sgbc_boundaries`)

The optional `sgbc_boundaries` object controls the boundary conditions at the two ends of the physical layer stack (`left` is the Elem1 / field-side interface; `right` is the Elem2 / inner or backing interface):

- **Passthrough sheet (interior face):** Omit `sgbc_boundaries` (or leave both sides coupled). Waves impinging from either side of the interior face reflect and transmit through the layer stack into the opposite volume.
- **PEC-backed coating (e.g. RAM on a metallic airframe):** Set `"sgbc_boundaries": { "right": "PEC" }`. The right side of the stack is terminated in a short circuit ($\mathbf{E}_t = 0$), and no field is transmitted past the coating.

---

## Automatic 1D layer meshing and order selection

You only need to specify the physical properties (`relative_permittivity`, `relative_permeability`, `bulk_conductivity`) and thickness (`material_width`) of each layer. At startup, `driver.cpp` inspects the excitation source's maximum frequency $f_{\max}$ and $1/e$ power frequency $f_{1e}$ to discretize each layer automatically:

1. **Wavelength, skin depth, and Péclet number:**
   For a layer with $\varepsilon = \varepsilon_r\varepsilon_0$, $\mu = \mu_r\mu_0$, and $\sigma = \sigma_{\mathrm{SI}}$, the solver computes:

$$
\lambda = \frac{1}{f_{\max}\sqrt{\mu\varepsilon}}, \qquad \delta = \frac{1}{\sqrt{\pi f_{\max}\mu\sigma}}, \qquad \mathrm{Pe} = \frac{\sigma}{f_{\max}\varepsilon}.
$$

2. **Adaptive polynomial order (`order`):**
   - High-loss / conductor regime ($\mathrm{Pe} > 100$): `order = 2` (captures steep exponential skin-depth decay without high-order ringing).
   - Low-loss dielectric regime ($\mathrm{Pe} < 0.1$ and $\sigma < 10^{-3}\,\mathrm{S/m}$): `order = 4`.
   - Moderate-loss regime: `order = 3`.

3. **Automatic segment count (`num_of_segments`):**
   Targets at least $20$ DOFs per wavelength ($\Delta x_{\mathrm{wave}} = \lambda / \lceil 20 / p \rceil$) and at least $2$ elements per skin depth ($\Delta x_{\mathrm{skin}} = \delta / 2$):

$$
N_{\mathrm{seg}} = \mathrm{clamp}\!\left(\left\lceil \frac{d}{\min(\Delta x_{\mathrm{wave}},\, \Delta x_{\mathrm{skin}})} \right\rceil,\, 4,\, 1000\right).
$$

Both `order` and `num_of_segments` can be overridden manually inside any layer object if desired. If a layer exceeds $7$ skin depths at $f_{1e}$ ($< -60\,\mathrm{dB}$ one-way transmission), the solver prints a startup notice suggesting a `"PEC"` boundary instead.

---

## Sub-stepping, `sgbc_cfl`, and LU caching

The number of SDIRK33 sub-steps per parent RK4 stage is chosen from the total acoustic crossing time $T_{\mathrm{cross}} = \sum_i d_i \sqrt{\mu_i\varepsilon_i}$ of the layer stack:

$$
\delta t_{\mathrm{rec}} = T_{\mathrm{cross}} \times \mathrm{sgbc\_cfl} \times \mathrm{opacity\_relax},
$$

where:
- `sgbc_cfl` defaults to `0.5` and can be set under `solver_options.sgbc_cfl` ($> 0$).
- For optically thick conductive layers where $d / \delta_{1e} > 3$ (evaluated at the source $1/e$ power frequency $f_{1e}$), $\mathrm{opacity\_relax} = \min\!\big((d/\delta_{1e})^2,\, 50\big)$ relaxes the sub-step size because diffusion dominates wave propagation.

Because RK4 uses stage increments of $0.5\Delta t$ (stages 1–3) and $1.0\Delta t$ (stage 4), each SGBC sub-solver encounters at most two distinct sub-step sizes $\delta t_{1/2}$ and $\delta t_1$. `GlobalEvolution::ImplicitSolve` caches the dense LU factorizations of $J = I - \gamma_{\mathrm{SDIRK}}\,\delta t\,A_{\mathrm{SGBC}}$ in a two-slot cache on the first time step and reuses them for the rest of the simulation without re-factorizing.

---

## JSON configuration

SGBC entries are defined in `model.boundaries` with `"type": "SGBC"` and either a single `"material"` object or a `"layers"` array (ordered from `left` to `right`).

| Field | Location | Default | Description |
|-------|----------|---------|-------------|
| `tags` | boundary object | — | Mesh boundary/interface attribute IDs. |
| `type` | boundary object | — | `"SGBC"`. |
| `material` | boundary object | — | Single-layer definition (mutually exclusive with `layers`). |
| `layers` | boundary object | — | Array of layer objects ordered from `left` (field side) to `right` (backing side). |
| `sgbc_boundaries` | boundary object | — | Optional `{ "left": "PEC"|"PMC"|"SMA", "right": "PEC"|"PMC"|"SMA" }`. |
| `relative_permittivity` | layer / material | `1.0` | Layer $\varepsilon_r$ ($\ge 1.0$). |
| `relative_permeability` | layer / material | `1.0` | Layer $\mu_r$ ($\ge 1.0$). |
| `bulk_conductivity` | layer / material | — | Required layer conductivity $\sigma_{\mathrm{SI}}$ in S/m ($\ge 0$). |
| `material_width` | layer / material | — | Required layer thickness $d$ in metres ($> 0$). |
| `num_of_segments` | layer / material | auto | Optional manual override for segment count in this layer. |
| `order` | layer / material | auto | Optional manual override for FE polynomial order in this layer. |
| `sgbc_cfl` | `solver_options` | `0.5` | Crossing-time CFL multiplier ($> 0$) for SGBC sub-stepping. |

### Example: PEC-backed 3-layer RAM coating

From [`testData/maxwellInputs/1D_SGBC_Multilayer/1D_SGBC_Multilayer.json`](../testData/maxwellInputs/1D_SGBC_Multilayer/1D_SGBC_Multilayer.json):

```json
{
  "tags": [3],
  "type": "SGBC",
  "sgbc_boundaries": {
    "right": "PEC"
  },
  "layers": [
    {
      "material_width": 0.0005,
      "relative_permittivity": 3.0,
      "relative_permeability": 1.0,
      "bulk_conductivity": 0.0
    },
    {
      "material_width": 0.003,
      "relative_permittivity": 12.0,
      "relative_permeability": 1.0,
      "bulk_conductivity": 2.0
    },
    {
      "material_width": 0.0005,
      "relative_permittivity": 4.0,
      "relative_permeability": 1.0,
      "bulk_conductivity": 10.0
    }
  ]
}
```

### Example: 3D coated sphere for RCS

From [`testData/maxwellInputs/3D_RCS_SGBC_CoatedSphere_O2/3D_RCS_SGBC_CoatedSphere_O2.json`](../testData/maxwellInputs/3D_RCS_SGBC_CoatedSphere_O2/3D_RCS_SGBC_CoatedSphere_O2.json):

```json
{
  "tags": [1],
  "type": "SGBC",
  "material": {
    "relative_permittivity": 4.0,
    "relative_permeability": 1.0,
    "bulk_conductivity": 0.05,
    "material_width": 0.05,
    "num_of_segments": 6,
    "order": 2
  },
  "sgbc_boundaries": {
    "right": "PEC"
  }
}
```

---

## ParaView export of internal SGBC fields

When top-level `probes.exporter` is enabled in the case JSON, `ProbesManager` automatically exports the internal 1D sub-solver fields for each SGBC boundary tag alongside the main volume ParaView collection:

```text
exports/ParaView/<run-mode>/<case>.msh/
├── <case>.msh.pvd
└── <case>_sgbc_tag<tag>.msh/
    └── <case>_sgbc_tag<tag>.msh.pvd
```

This allows you to inspect the standing-wave and skin-depth profiles inside the physical layers directly in ParaView.
