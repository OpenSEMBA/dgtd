# Session — PML sign audit and t = 20 stability (2026-05-28)

Follow-up to [10-session-2026-05-28-stability-dft.md](./10-session-2026-05-28-stability-dft.md). Plan: sign/discretization audit (do **not** cap `final_time` as the fix).

## Code added

| Item | Location |
|------|----------|
| `PML_SIGN_TEST` env (`0`–`8`) | [`PMLSignTest.h`](../../src/components/PMLSignTest.h), [`PMLSignTest.cpp`](../../src/components/PMLSignTest.cpp) |
| Sign modes S1–S6 | [`DGOperatorFactory.h`](../../src/components/DGOperatorFactory.h) `collectPMLOperatorBlocks` / `collectPMLComponentDriverBlocks` |
| S5 global mult sign | [`GlobalEvolution.cpp`](../../src/evolution/GlobalEvolution.cpp) `getPMLOperatorMultSign()` |
| Outer / inner PML `max\|psi\|`, `outer_Ey` | [`GlobalEvolution.cpp`](../../src/evolution/GlobalEvolution.cpp) `initPMLLocalizationDiagnostics`, `[PML Mult diag]` |
| PEC outer case | [`1D_PML_PEC/`](../../testData/maxwellInputs/1D_PML_PEC/) |
| dt-half reference case | [`1D_PML_dtHalf/`](../../testData/maxwellInputs/1D_PML_dtHalf/) (`time_step = 0.0005`) |
| Baseline outer BC | [`1D_PML.json`](../../testData/maxwellInputs/1D_PML/1D_PML.json) tag 2 → **SMA** |

### Sign-test matrix

| ID | `PML_SIGN_TEST` | Change | t = 20 `max\|psi\|` (dt = 0.001) |
|----|-----------------|--------|-----------------------------------|
| default | 0 | Current conventions | **1.77×10¹³** |
| S1 | 1 | Flip H/E correction signs | 6.30×10¹⁰⁵ (worse) |
| S2 | 2 | Face driver +w (no SBP) | 2.43×10²⁷ (worse) |
| S3 | 3 | Face cross-column (global layout) | 7.00×10⁴⁷ (worse) |
| S4 | 4 | Negate driver weight w | 6.30×10¹⁰⁵ (same as S1) |
| S5 | 5 | `PMLOperator_->AddMult(..., -1)` | **1.77×10¹³** (same as default) |
| S6 | 6 | ψ mass +1 | **1.77×10¹³** (α = 0 → mass inactive) |
| S7 | 7 | Terminating boundary face on ψ driver | **2.18×10¹⁵** (worse; opt-in only) |

**Criterion:** No test achieved bounded `max\|psi\|` (O(1–10²)) at t = 20.

## Localization (default, dt = 0.001, SMA @ x = 3)

| Approx. time | Mult call | max \|ψ\| | outer PML | max \|Ey\| | outer Ey |
|--------------|-----------|------------|-----------|------------|----------|
| 7.5 s | 30000 | 2.4 | — | ~1 | — |
| 11.25 s | 45000 | 0.17 | — | ~1 | — |
| **13.75 s** | **55000** | **1417** | onset | ~1 | — |
| 12 s | 48000 | **2.43** | outer ≈ global | ~0.08 | — |
| 20 s | 80000 | **1.77×10¹³** | **1.48×10¹³** | **1.15×10¹²** | **1.15×10¹²** |

ψ-led growth starts in **outermost PML element** after the bulk pulse has passed; fields eventually diverge with ψ.

## PEC vs SMA at x = 3 (tag 2)

| Outer BC | t = 20 max \|ψ\| |
|----------|------------------|
| SMA (baseline) | 1.77×10¹³ |
| PEC ([`1D_PML_PEC`](../../testData/maxwellInputs/1D_PML_PEC/)) | **4.54×10²⁷** |

Blow-up **persists and worsens** with PEC → not SMA reflection alone; primary issue is PML operator / explicit stiff coupling, not only global BC mismatch.

## dt refinement (not a sign fix)

[`1D_PML_dtHalf`](../../testData/maxwellInputs/1D_PML_dtHalf/) (`dt = 0.0005`, t = 20):

| Mult call | Time | max \|ψ\| |
|-----------|------|------------|
| 80000 | 10 s | 0.132 |
| 160000 | 20 s | **1.77×10¹³** |

Halving `dt` **delays** onset (~10 s stable) but **does not** fix t = 20 under explicit RK4.

## Split audit (sign sweep failed)

### Gedney ↔ CSR map (1D TE, stretch X)

| Continuous | Discrete block | Sign (default) |
|------------|----------------|----------------|
| \(\dot\psi^E_{X,Z} = -\alpha\psi + \sigma\mathcal{D}_x(E_y)\) | ψ mass + driver | mass **−1**, vol **+w**, face **−w** on Ey col |
| \(\dot H_z \mathrel{-}= \psi^E_{X,Z}/\kappa\) | H[z] ← ψ^E | **−1** |
| \(\dot\psi^H_{X,Y} = -\alpha\psi + \sigma\mathcal{D}_x(H_z)\) | ψ^H driver | same pattern on Hz col |
| \(\dot E_y \mathrel{+}= \psi^H_{X,Y}/\kappa\) | E[y] ← ψ^H | **+1** |
| Global \(\dot H_z\) from \(\partial_x E_y\) | `globalOperator_` | full domain including PML |

### Duplicate-derivative hypothesis

Architecture ([`02-codebase-architecture.md`](./02-codebase-architecture.md)): `globalOperator_` applies the **full** curl in PML; `PMLOperator_` **adds** σ-weighted \(\mathcal{D}_d\) into ψ and ψ/kappa into fields. With **κ = 1**, the κ-part of stretch is in the global curl; ψ should carry the CFS memory. A wrong **relative** sign between global face flux and PML face SBP at the **terminating** layer (interior-only ψ face integrator) remains a suspect, but trial **boundary face** terms (S7) **increased** growth.

**Next implementation candidates (not done here):**

1. Mark PML elements in `collectGlobalOneNormalOperators` / directional blocks and subtract the stretch-direction contribution from global curl where ψ driver duplicates it.
2. Milestone B: implicit RK / SDIRK with extended `ImplicitSolve` for stiff ψ modes.
3. Revisit cross-column face placement with correct row/column `FieldOffsets` (S3 trial used alt-field **columns** only; may need full `(f[y], alt[f][z])` row–col pairs).

## Regression

| Case | Result |
|------|--------|
| `1D_PEC` | OK (no extended state) |
| `1D_PML` t ≈ 12 s, dt = 0.001 | max \|ψ\| ≈ **2.4** (acceptable) |
| `1D_PML` t = 20 s, dt = 0.001 | **Not stable** |

## Conclusion

**No single sign flip (S1–S6, S5 global mult) fixes late-time blow-up.** Diagnostics confirm **outer-PML ψ** growth inward after ~t ≈ 13–14. Stability to t = 20 likely needs **discretization split** (global vs PML derivative partitioning) and/or **implicit/stiff time integration**, not a TFSF-style global minus on `PMLOperator_`.

Sign conventions in [`01-physics-and-formulation.md`](./01-physics-and-formulation.md) remain **unchanged** (defaults retained).
