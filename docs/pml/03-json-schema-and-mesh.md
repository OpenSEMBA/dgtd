# JSON schema, Gmsh mesh, and stretch depth

## Design principle

**The mesh defines geometry and thickness. JSON defines material type and grading law.**

There is **no required `pml_thickness` field** in JSON. The user models a PML **volume in Gmsh**, assigns an **element attribute tag**, and references that tag in JSON with `"type": "PML"`.

---

## Workflow (user-facing)

```text
1. Gmsh: draw vacuum + PML volumes
2. Gmsh: assign Physical Volume tags (attributes)
3. Export .msh
4. JSON: map tags to "vacuum" or "PML" with stretch parameters
5. dgtd: at init, compute profiles from mesh geometry + JSON
6. Run with probes in vacuum; DFT offline for −40 dB check
```

---

## Material types in JSON

### Vacuum

```json
{
  "tags": [9, 10],
  "type": "vacuum"
}
```

Semantics:

- ε = 1, μ = 1, σ_bulk = 0 (normalized)
- Equivalent to omitting material keys in legacy JSON

**Implementation:** map to `buildVacuumMaterial()` in `GeomTagToMaterial`.

### PML (volumetric CFS-CPML)

```json
{
  "tags": [1, 2],
  "type": "PML",
  "matches_vacuum": true,
  "grading_order": 4,
  "target_reflection": 1e-6,
  "kappa_max": 1.0,
  "alpha_max": 0.0,
  "active_axes": ["X"]
}
```

| Field | Type | Required | Semantics |
|-------|------|----------|-----------|
| `tags` | int[] | yes | Gmsh volume (or region) attribute IDs |
| `type` | string | yes | Must be `"PML"` |
| `matches_vacuum` | bool | yes (default true) | ε = μ = 1 in PML; stretch provides absorption only |
| `grading_order` | int | yes | Power-law exponent **m** for profile vs normalized depth |
| `target_reflection` | double | yes | Design reflection level (e.g. 1e-6); used to set **σ_max** (implementation formula TBD from Taflove/Gedney) |
| `kappa_max` | double | yes | Maximum κ at outer PML edge (≥ 1); 1.0 = no κ grading in tests |
| `alpha_max` | double | yes | Maximum α for CFS; **0** = classical PML limit; use **> 0** for late-time stability in production |
| `active_axes` | string[] | yes | Subset of `"X"`, `"Y"`, `"Z"`; which directions are stretched for **this tag block** |

**Forbidden on PML tags:**

- `bulk_conductivity`
- `relative_permittivity` / `relative_permeability` unless explicitly allowed later (currently **matches_vacuum** only)

**Multiple PML blocks:** Reference case `2D_RCS_Circle_Vol_PML.json` uses separate blocks for X-only, Y-only, and XY corner tags — same pattern user should follow.

---

## Reference case structure

Path: `testData/maxwellInputs/2D_RCS_Circle_Vol_PML/2D_RCS_Circle_Vol_PML.json`

| Tag block | active_axes | Role |
|-----------|-------------|------|
| 9, 10 | (vacuum) | Interior / scatterer region |
| 1, 2 | X | X-directed PML slabs |
| 3, 4 | Y | Y-directed PML slabs |
| 5–8 | X, Y | Corner PML |

User will add **1D acceptance case** under `testData/maxwellInputs/1D_PML/`.

Suggested 1D JSON skeleton:

```json
{
  "solver_options": {
    "order": 2,
    "time_step": 0.01,
    "final_time": 4.0,
    "upwind_alpha": 1.0
  },
  "model": {
    "filename": "1D_PML_Slab.msh",
    "materials": [
      { "tags": [1], "type": "vacuum" },
      {
        "tags": [2],
        "type": "PML",
        "matches_vacuum": true,
        "grading_order": 4,
        "target_reflection": 1e-4,
        "kappa_max": 1.0,
        "alpha_max": 0.0,
        "active_axes": ["X"]
      }
    ],
    "boundaries": [
      { "tags": [1], "type": "SMA" },
      { "tags": [2], "type": "SMA" }
    ]
  },
  "sources": [ ... ],
  "probes": { ... }
}
```

Adjust tags, sources, boundaries to match user's Gmsh file. Outer boundaries: **SMA** preferred over PEC for absorber tests (see locked decisions).

### `upwind_alpha` and PML

| Value | Effect |
|-------|--------|
| **`0.0`** | Centered flux only (OneNormal + volume). **Recommended** for late-time ψ stability in current 1D matrix. |
| **`1.0`** | Full Hesthaven upwind (Zero/Two scaled by `upwind_alpha`). `PMLOperator_` ψ driver includes matching σ-weighted Zero/Two on PML interior faces when α > 0 — see [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md). |

Reference cases: [`1D_PML_centered/`](../../testData/maxwellInputs/1D_PML_centered/), [`1D_PML_buffer_centered/`](../../testData/maxwellInputs/1D_PML_buffer_centered/). Baseline upwind: [`1D_PML/`](../../testData/maxwellInputs/1D_PML/), [`1D_PML_buffer/`](../../testData/maxwellInputs/1D_PML_buffer/).

Uniform probes (1.99 / 2.01 / 2.99) and stability gate: [`check_pml_probe_stability.py`](../../testData/maxwellInputs/1D_PML/check_pml_probe_stability.py).

**Note:** Runtime **`PML_DERIV_SPLIT`** (marker subtraction from `globalOperator_`) was **removed** — it broke the vacuum–PML interface. Details: [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md).

---

## How thickness works without JSON

Example: 1D mesh, vacuum tag 1 on [0, 0.8], PML tag 2 on [0.8, 1.0].

- Physical PML thickness = **0.2** (mesh coordinates, normalized)
- At QP with x = 0.9 in PML element: depth **ρ_X = 0.1** from interface at x = 0.8
- At QP with x = 1.0: **ρ_X = 0.2** = local layer depth for that element

**No JSON number duplicates 0.2.**

---

## Computing depth ρ_d at quadrature points

### Goal

For each PML element **K**, each quadrature point **x_q**, each active direction **d**:

- **ρ_d(x_q)** = distance from **x_q** into the PML along direction **d**, measured from the **vacuum–PML interface** normal to **d**.

### Algorithm (Milestone A — attribute neighbor walk)

**Init phase (once):**

1. Build set **PMLTags** from all `PMLProperties.geom_tags`.
2. For each element `el`, `attr = mesh.GetAttribute(el)`; `is_pml = attr ∈ PMLTags`.
3. For each PML element `el` and each QP `x_q`:
   - For each active direction **d** (from properties for `attr`):
     - Find **interface depth**: minimum distance from `x_q` to points on faces shared with a **non-PML** neighbor, measured along **d** (sign: positive into PML).
     - **Practical serial approach:**
       - Use element centroid and bounding box along **d** relative to mesh-wide vacuum/PML partition.
       - Or ray-walk: from `x_q`, step toward **+d** until leaving PML region; interface is where attribute changes.

**Corner tags (active X and Y):**

- Compute **ρ_X** and **ρ_Y** **independently** using the same rule per direction.
- Do **not** require diagonal distance.

### Interface matching (profiles at ρ = 0)

Enforce:

```text
κ_d(0) = 1,  σ_d(0) = 0,  α_d(0) = 0
```

Power-law toward max at outer edge of tagged region along **d**:

```text
ξ = ρ_d / L_d   ∈ [0, 1]
σ_d(ρ) = σ_d,max * ξ^m
κ_d(ρ) = 1 + (κ_max - 1) * ξ^m     (if κ_max > 1)
α_d(ρ) = α_max * ξ^m
```

where **m = grading_order**, **L_d** = max depth in PML along **d** for that element (or local distance to outer boundary of PML tag region).

### σ_max from target_reflection

**Open implementation task:** use Taflove CPML or Gedney guidance to relate single-reflection estimate to **σ_d,max** given **L_d**, **m**, and `target_reflection`.

Placeholder policy until tuned:

```text
σ_d,max = -(m + 1) * ln(target_reflection) / (2 * L_d)
```

(normalized units; verify against literature and −40 dB acceptance.)

---

## Vacuum–PML interface (mesh + DG)

- Interface is an **interior face** between vacuum attribute and PML attribute.
- **No special boundary condition tag** required on that face.
- DG **interior flux** connects vacuum and PML elements.
- Stretch zeros at interface ensure continuous impedance matching in continuous limit.

**Do not** mark vacuum–PML interface as PEC/SMA/SGBC.

---

## Outer boundary of PML shell

Where volumetric PML meets the **exterior** mesh boundary (the terminating face of the outermost PML element):

| JSON `type` | Effect in `globalOperator_` |
|-------------|------------------------------|
| **`PML_NONE`** | No `AddBdrFaceIntegrator` on that tag — PML layer terminates without SMA/PEC flux (use when PML ends on ∂Ω). |
| **`PEC`** | Gedney-style PEC-backed slab (literature default for FDTD/DGFETD validation). |
| **`SMA`** | Standard dgtd absorbing BC — **avoid on the same face as terminating PML** (conflicts with PML ADE on shared DOFs). |

Reference cases:

- [`1D_PML/`](../../testData/maxwellInputs/1D_PML/) — SMA on outer tag (baseline A).
- [`1D_PML_PML_NONE/`](../../testData/maxwellInputs/1D_PML_PML_NONE/) — `PML_NONE` on outer tag (case B).
- [`1D_PML_PEC/`](../../testData/maxwellInputs/1D_PML_PEC/) — PEC on outer tag (case C).
- [`1D_PML_buffer/`](../../testData/maxwellInputs/1D_PML_buffer/) — vacuum buffer + SMA past PML (case D).

**Do not** mark the vacuum–PML **interface** (interior face) as PEC/SMA/SGBC/`PML_NONE`.

---

## Parser changes (`driver.cpp`)

In `assembleAttributeToMaterial()`:

1. Read `type` field per material block.
2. If `"vacuum"` → `buildVacuumMaterial()` for listed tags.
3. If `"PML"` → validate fields; append `PMLProperties`; **do not** insert into `GeomTagToMaterial` as conductive material.
4. If legacy block (no type) → keep current eps/mu/sigma parsing.

Validation errors:

- PML tag overlaps vacuum definition
- `active_axes` contains direction ≥ mesh dimension
- Unknown axis string
- `bulk_conductivity` on PML block

---

## README / schema documentation

Update `README.md` materials section:

- Document `"type": "PML"` fields
- Remove `SBC_PML` boundary documentation
- Cross-link `docs/pml/`

Note: `.cursor/rules` mention `core/schema` — if JSON schema is added later, mirror this document.

---

## FAQ: confusion resolved in planning

**Q: If mesh defines thickness, why JSON grading params?**  
A: Mesh gives **L**; JSON gives **how fast** κ, σ, α grow over **L** and their maxima / design reflection.

**Q: Could we add optional `pml_thickness`?**  
A: Redundant with Gmsh; skip for v1 unless auto-meshing without CAD requires it later.

**Q: Does Paraview need PML tags?**  
A: Exporter shows mesh attributes as today; ψ is internal only.
