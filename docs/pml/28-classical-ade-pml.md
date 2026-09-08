# Classical ADE-PML (CuDG3D-style) — session 28 (2026-09-07)

**Status: active path** after Gedney CFS pause ([`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md)).

## Formulation

Uniaxial or multi-axis stretch (normalized \(\varepsilon=\mu=1\)), volume-only. For each active stretch axis \(s\) on a material block (historical CuDG3D ADE, \(\alpha=0\)):

```text
∂t E_s  +=  σ_s E_s − J_s
∂t E_⊥  −=  σ_s E_⊥
∂t J_s   =  σ_s² E_s − σ_s J_s
(same for H / M)
```

The stretch-parallel \((E_s,J)\) block is only **marginally** stable and can leave late-time wall modes (Ex in X-PML, Ey in Y-PML). CFS-lite `alpha_max` is **not** part of this path (rejected in JSON with `kappa_max`).

Biaxial/triaxial ADE is the **superposition** of those stacks when `active_axes` lists more than one direction. No face ADE (ordinary DG fluxes at vacuum–PML). SMA on the outer PML face does **not** remove volumetric \(\parallel\)-ADE modes.

**σ discretization:** `MassIntegrator` samples σ(x) (and σ²(x)) at each quadrature point through `PMLProfileCoefficient` → [`PMLProfileData::evaluateAtTransform`](../../src/components/PMLProfiles.h) (power-law in depth for `grading_order >= 1`; constant for `0`). Element-mean flattening was removed — it mismatched graded σ² and could distort H relative to E.

**Depth model:**
- **`stretch_mode: "box"`** (default): \(\rho_s\) / \(L_s\) from planar vacuum–PML interfaces per axis. Both **+** and **−** sides of each axis are tracked (e.g. top and bottom Y-slabs); previously only one sign was kept, which zeroed graded \(\sigma\) on the opposite slab.
- **`stretch_mode: "radial"`**: \(\rho=\max(0,\|x-c\|-r_{\mathrm{in}})\), \(L=r_{\mathrm{out}}-r_{\mathrm{in}}\); same \(\sigma(\rho)\) feeds every listed `active_axes` stack (**radial σ grading, Cartesian ADE** — not full cylindrical UPML). Optional `radial_center`; otherwise inferred from interface face centers.

**Blocks:** one tag with `["X","Y"]` is supported; separate uniaxial blocks (as in `2D_RCS_Circle_Vol_PML`) remain valid.
## State and Mult

```text
[Ex Ey Ez Hx Hy Hz | J_d0 M_d0 J_d1 M_d1 … ]
n_aux = 2 × (# stretch dirs) × ndofs
```

Layout: [`ClassicalPMLLayout`](../../src/components/ClassicalPMLLayout.h).

`Mult()` order:

```text
SGBC → pack EH(+nbr) → globalOperator_ → classicalPMLOperator_ AddMult → TFSF → SGBC flux
```

Assembled operator: `DGOperatorFactory::buildClassicalPMLOperator` — \(M^{-1}M(\sigma)\) field blocks + \(J\)/\(M\) couplings on PML markers.

**MPI:** Profiles are built on the full serial mesh (global vacuum–PML interface and \(L\)). `evaluateAtTransform` keys off `T.Attribute` + physical coordinates — not parallel `ElementNo`. ADE is volume-local (no aux ghost exchange). Each rank owns a local ODE of size \(6N_{\mathrm{loc}}+n_{\mathrm{aux}}\). Ranks with no local PML elements omit the classical operator (`nullptr` no-op in `Mult`).

## JSON (unchanged region contract)

```json
{
  "tags": [3],
  "type": "PML",
  "matches_vacuum": true,
  "grading_order": 4,
  "target_reflection": 1e-6,
  "active_axes": ["X", "Y"],
  "stretch_mode": "box"
}
```

| Field | Default | Meaning |
|-------|---------|---------|
| `stretch_mode` | `"box"` (`0`) | `"box"` / `0`: planar depth per axis. `"radial"` / `1`: \(\rho=\max(0,\|x-c\|-r_{\mathrm{in}})\); ADE stacks stay Cartesian (radial σ profile only). |
| `radial_center` | inferred | Optional `[x,y]` / `[x,y,z]`. Only with `stretch_mode: "radial"`. If omitted, mean of vacuum–PML interface face centers. |

**Rejected:** `kappa_max`, `alpha_max` (CFS-only).

## Acceptance

| Case | Gate |
|------|------|
| `1D_PML` | DFT reflection ≤ −40 dB ([`scripts/pml_dft_reflection.py`](../../scripts/pml_dft_reflection.py)) |
| `2D_PML_X_slab` | Stable finite fields through `final_time` |
| `2D_RCS_Circle_Vol_PML` | Assembles (X+Y uniaxial blocks) |
| `2D_RCS_Circle_1m_G2_PML` | Assembles biaxial onion-ring with `stretch_mode: "radial"` |
