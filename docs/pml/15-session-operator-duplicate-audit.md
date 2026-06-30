# Session — Operator assembly and duplicate-curl audit (2026-05-27)

Follow-up to [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md). Implements the attached **PML operator duplicate audit** plan (diagnosis only).

## Runtime flags

| Env | Effect |
|-----|--------|
| `PML_OPERATOR_AUDIT=1` | Assembly manifest, Frobenius overlap, CSR export under `Exports/Operators/<case>_audit/`, integrator baseline print |
| `PML_MULT_PROBE=1` | After operator build: isolated `globalOperator_->Mult` + `applyPMLCoupling` probes (no TFSF/SGBC) |
| `PML_SKIP_GLOBAL=1` | Skip `globalOperator_->Mult` in `Mult()` (audit only) |
| `PML_SKIP_PML=1` | Skip `PMLOperator_->AddMult` in `Mult()` (audit only) |

Code: [`PMLOperatorAudit.h`](../../src/components/PMLOperatorAudit.h), [`PMLOperatorAudit.cpp`](../../src/components/PMLOperatorAudit.cpp), hooks in [`DGOperatorFactory.h`](../../src/components/DGOperatorFactory.h), [`GlobalEvolution.cpp`](../../src/evolution/GlobalEvolution.cpp).

---

## Phase 1 — Integrator baseline (1D TE, `active_axes: ["X"]`)

| Path | Integrator | Domain | Coefficient | Maps into (1D TE) |
|------|------------|--------|-------------|---------------------|
| Global | `DerivativeIntegrator` + `MInv[H]` | Full mesh | 1 | `Ħ_z ← ∂_x E_y` (vol, weight −1) |
| Global | `MaxwellDGOneNormalJumpIntegrator` + `MInv[H]` | All interior + BCs | 1 | Same `H[z]←E[y]` SBP face flux (weight −1) |
| Global | `MaxwellDGZeroNormal` / `TwoNormal` | Full mesh | `upwind_alpha` | Self-blocks; **inactive when α=0** |
| PML ψ driver | `DerivativeIntegrator(σ)` + scalar `MInv` | PML marker vol | σ(x) | `ψ̇^E_{X,Z} ← σ ∂_x E_y` (vol +w) |
| PML ψ driver | `MaxwellDGCoefficientOneNormalJumpIntegrator(σ)` | PML interior faces | σ(x) | SBP face on `E[y]` (−w on same column) |
| PML ψ driver | σ×α Zero/Two | PML marker | σ·α | Only when `upwind_alpha > 0` |
| PML correction | `MassIntegrator(1/κ)` + `MInv[H]` | PML vol | 1/κ | `Ħ_z ← −ψ^E/κ` (default sign) |

`pmlCurlDerivWeight(H, Z, E, Y, X) = +1` matches global directional placement for the active 1D chain.

---

## Phase 2 — Assembly audit (`1D_PML` vs `1D_PML_centered`)

Command (from case directory):

```bash
PML_OPERATOR_AUDIT=1 PML_MULT_PROBE=1 \
  ../../../../build/gnu-release-mpi/bin/opensemba_dgtd \
  -i $(pwd)/1D_PML.json
```

### Global operator manifest

| Metric | `1D_PML` (α=1) | `1D_PML_centered` (α=0) |
|--------|----------------|-------------------------|
| Block placements after Directional | 4 | 4 |
| After OneNormal | 8 | 8 |
| After ZeroNormal | 14 | 14 (integrator coeff 0 → **no effective nnz**) |
| After TwoNormal | 16 | 16 (same) |
| **Merged global nnz** | **6192** | **3360** |
| `||Hz←Ey||_F` (full global) | 776.46 | 776.46 |
| `||Hz←Ey||_F` on PML element DOFs only | 314.55 | 314.55 |
| `global_Dx+OneNormal_x` `||Hz←Ey||_F` | 1069.06 | 1069.06 |
| Ratio Dx+OneNormal / full | 1.38 | 1.38 |

**Interpretation:** The Ey→Hz semidiscrete block is already dominated by **directional + OneNormal** in x (not Zero/Two). Centered vs upwind differs only in **extra global nnz** from Zero/Two placements (6192 vs 3360); the Ey→Hz block norm is **unchanged** at α=0.

### PML operator manifest

| Metric | `1D_PML` (α=1) | `1D_PML_centered` (α=0) |
|--------|----------------|-------------------------|
| PMLOperator nnz | 828 | 936 |
| `driver_upwind_nnz` | 432 | 0 |
| `driver_vol_nnz` / `driver_face_nnz` | 360 / 432 | same |
| `||ψ←Ey||_F` (assembled driver, vol−face SBP) | **3677.07** | **3677.07** |
| `||Hz←ψ||_F` (ψ→H correction block) | 5.0 | 5.0 |
| `||ψ←ψ_mass||_F` (α=0 in JSON) | 0 | 0 |

**Interpretation:** The σ-weighted ψ driver is **much larger** in Frobenius norm than the global `Hz←Ey` restriction to PML DOFs (314.55), even though both are built from the same weak `D_x` + face SBP pattern. Contributors:

1. **σ(x) multiplies** every PML quadrature point (large in the layer).
2. **Scalar** `MInv` on ψ rows vs **Maxwell** `MInv[H]` on global curl.
3. Driver is a **different state row** (ψ) not directly comparable to `Hz←Ey` without the ψ ODE chain.

Exported CSRs: `Exports/Operators/1D_PML_audit/` and `1D_PML_centered_audit/` (`global_full.csr`, `global_Dx_only.csr`, `pml_driver_Ey_to_psi.csr`, `pml_full.csr`).

---

## Phase 3 — Mult() probes

| Probe | Input | `max|dHz/dt|` global | after PML | `max|dpsi/dt|` |
|-------|-------|----------------------|-----------|----------------|
| P1 vacuum `E_y` | unit on one vacuum Ey DOF | 0 | 0 | 0 |
| P2 PML `E_y` | unit on one inner PML Ey DOF | 0 | 0 | **0.0032** |
| P3 unit `ψ^E_{X,Z}` | unit on one ψ DOF | 0 | 0 | 0 |
| P4 zero | — | 0 | 0 | 0 |

**Notes:**

- Probes call **`globalOperator_->Mult` + `applyPMLCoupling` only** (no TFSF/SGBC), with a single local DOF excited and **no face-neighbor exchange** — DG flux therefore does not produce `Ħ_z` from an isolated `E_y` pulse (P1/P2 `Hz` = 0 is expected for this discrete test, not a wiring bug).
- **P2** confirms the **ψ driver**: nonzero `dpsi/dt` from unit `E_y` in PML.
- **P3** `Hz` from unit ψ may require a distributed ψ mode or RK mass term; static single-DOF test is inconclusive for correction magnitude.

`multWorkVec` vs `pmlWorkVec` field blocks: **max|diff| = 0** at first `applyPMLCoupling`.

---

## Conclusions (ranked)

| # | Finding | Confidence |
|---|---------|------------|
| 1 | **Duplicate-stretch hypothesis confirmed at architecture level:** `globalOperator_` applies full `D_x` + OneNormal curl in PML; `PMLOperator_` adds a **second** σ-weighted `D_x` + face chain into ψ, then ψ/kappa feeds back into `H`. This is **not** the Gedney split \(1/s = (1/\kappa)\partial + \text{memory}\) with memory-only in ψ. | High |
| 2 | **Centered vs upwind** does not change Ey→Hz block norms; upwind only adds global/PML Zero/Two nnz when α>0. Late-time centered stability (session 14) is **not** explained by matching this block alone. | High |
| 3 | **Mass model mismatch:** ψ driver uses **scalar** L2 `MInv(1)`; global curl uses **Maxwell** `MInv[H]` (ε/μ). Align if pursuing a single discrete symbol for \(\mathcal{D}_x\). | Medium |
| 4 | **Interface faces:** PML face integrators use `pml_marker` on interior faces only; vacuum–PML interface still has **global** face flux. No separate ψ face term at the interface (S7 trial worsened stability). | Medium |
| 5 | Single-DOF Mult probes **do not** validate `Hz` coupling strength; use distributed tests or full simulation probes. | High |

---

## Recommended next implementation (out of scope here)

1. **Regional assembly** was tried in session 16 and **reverted** — omitting global curl/flux in PML broke DG consistency; PML auxiliaries should remain **additive** on the full global operator.
2. Investigate duplicate-curl via **signs**, ψ–field coupling, or formulation split without removing global flux.
3. Optional \(1/\kappa\) global mass when `kappa_max > 1`, ψ-driver `MInv` alignment.

---

## Test plan checklist

- [x] Assembly report for `1D_PML` and `1D_PML_centered`
- [x] Frobenius Ey→Hz / ψ←Ey metrics
- [x] Mult probes P1–P4 (with documented DG limitation)
- [x] Written conclusions in this doc
