# Classical ADE-PML (CuDG3D-style) — session 28 (2026-09-07)

**Status: active path** after Gedney CFS pause ([`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md)).

## Formulation

Uniaxial stretch axis \(s\) (normalized \(\varepsilon=\mu=1\)), volume-only:

```text
∂t E_s  +=  σ_s E_s − J_s
∂t E_⊥  −=  σ_s E_⊥
∂t J_s   =  σ_s² E_s − σ_s J_s
(same for H / M)
```

No face ADE (ordinary DG fluxes at vacuum–PML). Source: commented CuDG3D [`PMLUniaxial.hpp`](../../external/Cudg3d/src/core/dg/dispersive/PMLUniaxial.hpp).

**σ discretization:** `MassIntegrator` samples σ(x) (and σ²(x)) at each quadrature point through `PMLProfileCoefficient` → [`PMLProfileData::evaluateAtTransform`](../../src/components/PMLProfiles.h) (power-law in depth for `grading_order >= 1`; constant for `0`). Element-mean flattening was removed — it mismatched graded σ² and could distort H relative to E.

**v1 scope:** one `active_axes` entry per PML material block. Multi-axis single tags error out — use separate uniaxial blocks (as in `2D_RCS_Circle_Vol_PML`).

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
  "active_axes": ["X"]
}
```

## Acceptance

| Case | Gate |
|------|------|
| `1D_PML` | DFT reflection ≤ −40 dB ([`scripts/pml_dft_reflection.py`](../../scripts/pml_dft_reflection.py)) |
| `2D_PML_X_slab` | Stable finite fields through `final_time` |
| `2D_RCS_Circle_Vol_PML` | Assembles (X+Y uniaxial blocks) |
