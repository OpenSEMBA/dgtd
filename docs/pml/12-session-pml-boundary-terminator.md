# Session — PML outer boundary terminator experiments (2026-05-28)

Follow-up to [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md). Tests the **SMA + terminating PML** hypothesis with a new **`PML_NONE`** boundary type and a **vacuum buffer** mesh.

## Hypothesis

`globalOperator_` applies **SMA** face flux on ∂Ω while `PMLOperator_` drives ψ and corrections on the **same outermost PML element** — competing termination models (not literal double SMA in `PMLOperator_`).

## Implementation

| Item | Location |
|------|----------|
| `BdrCond::PML_NONE` | [`Types.h`](../../src/components/Types.h) |
| Markers `pmlNoneMarker_` / `intpmlNoneMarker_` | [`Model.h`](../../src/components/Model.h), [`Model.cpp`](../../src/components/Model.cpp) |
| JSON `"PML_NONE"` | [`driver.cpp`](../../src/driver/driver.cpp) `assignBdrCond` |
| Skip `AddBdrFaceIntegrator` for `PML_NONE` | [`DGOperatorFactory.h`](../../src/components/DGOperatorFactory.h) `buildZeroNormal` / `buildOneNormal` / `buildTwoNormal` boundary loops |
| Case B | [`1D_PML_PML_NONE/`](../../testData/maxwellInputs/1D_PML_PML_NONE/) — tag 2 → `PML_NONE` |
| Case D | [`1D_PML_buffer/`](../../testData/maxwellInputs/1D_PML_buffer/) — +2 vacuum elements (x = 3.0–3.2), SMA at x = 3.2 |
| Schema | [`03-json-schema-and-mesh.md`](./03-json-schema-and-mesh.md) outer-boundary table |

**Note:** MFEM may emit `Non-positive attributes on the boundary!` when `PML_NONE` is used (no boundary bilinear form on that tag). Expected for this experiment.

## Comparison matrix (`final_time = 20`, `dt = 0.001`, `alpha_max = 0`)

| Case | Path | Outer BC | max \|ψ\| @ t≈7.5 | @ t≈12 | @ t=20 |
|------|------|----------|-------------------|--------|--------|
| **A** | `1D_PML/` | SMA @ x=3 | 2.40 | 2.43 | **1.77×10¹³** |
| **B** | `1D_PML_PML_NONE/` | PML_NONE @ x=3 | 2.40 | 2.43 | **1.77×10¹³** |
| **C** | `1D_PML_PEC/` | PEC @ x=3 | 4.78 | 3.25×10⁸ | **4.54×10²⁷** |
| **D** | `1D_PML_buffer/` | SMA @ x=3.2 (past PML) | 1.14 | **77** | **1.25×10¹⁵** |

Diagnostics: `[PML Mult diag]` calls 30000 / 48000 / 80000 (4 RK stages per step).

### Observations

1. **A vs B — identical ψ evolution** (diff only at float noise). Removing SMA boundary face integrators on the terminating PML face **does not** fix late-time blow-up or change mid-time ψ at t ≈ 12. The instability is **not explained by SMA face flux alone** on this 1D case.

2. **C (PEC)** — mid-time already degraded (ψ ~ 10⁸ at t ≈ 12); worst at t = 20. PEC on the same face as terminating PML is **not** the Gedney “PEC-backed slab” geometry (PEC is on the PML element face, not behind a separate buffer).

3. **D (vacuum buffer)** — PML still ends at x = 3 (interior); SMA at x = 3.2 on vacuum only. **Does not stabilize** t = 20 (still 10¹⁵). Mid-time **worse** than A/B at t ≈ 12 (max \|ψ\| ≈ 77 vs 2.4). Outer-element ψ at blow-up: `outer` ≈ 1.8×10¹⁰ ≪ global max — growth is **not confined** to the old outer PML cell only.

4. **Regression:** `1D_PEC` completes with no extended state.

## Conclusions

| Claim | Verdict |
|-------|---------|
| “Double SMA” in `PMLOperator_` | **Rejected** — PML does not add SMA integrators |
| SMA on ∂Ω conflicts with terminating PML | **Weak effect** on this case — `PML_NONE` ≈ SMA for ψ blow-up |
| Vacuum buffer + distant SMA fixes t = 20 | **Rejected** |
| PEC outer tag fixes t = 20 | **Rejected** (worse) |

**Recommended next work** (unchanged from sign audit):

1. **Derivative split** — remove or subtract stretch-direction curl in PML from `globalOperator_` where ψ driver duplicates \(\mathcal{D}_d\).
2. **Gedney §V weak-form** — terminating face / \(\hat{n}=\hat{\nu}/\kappa_\nu\) (not another SMA toggle).
3. **Time integration** — α > 0, smaller `dt`, or implicit ψ subsystem.

## Product guidance (interim)

- Use **`PML_NONE`** on outer tags when you want PML to be the only boundary model on that face (avoids SMA flux; MFEM may warn).
- Do **not** assume `PML_NONE` alone stabilizes late time on `1D_PML`.
- For literature-aligned termination, prefer **PEC behind a PML slab** with mesh geometry that places PEC **outside** the σ-profile region (not tested successfully here with tag-2 PEC on the PML element face).
- **Vacuum buffer** is a useful diagnostic mesh pattern but **not** a fix by itself.

## Visual correlation (user screenshots)

Transit-time peak sharpening at x = 3 (t ≈ 7.9 s) occurs for **A and B alike** — consistent with PML σ-profile + interior DG flux, not only SMA. Late border fuzz (t ≳ 11.5 s) matches ψ-led outer growth in diagnostics for A/B/C/D.
