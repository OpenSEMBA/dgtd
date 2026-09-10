# Gedney CFS-CPML attempt — paused

**Status: paused / not an active code path.**

## 2026-09-09 — second attempt removed

Reintroduced `pml_formulation: "gedney"` with shared centered \(D\) (Derivative+OneNormal) driving \(\psi\). Centered `1D_PML_Gedney` DFT failed (~+28 dB; late \(|E_y|\sim 1\)). Code path deleted again in favor of Bagci/Chen **SC-PML** field-driven ADE ([`30-sc-pml-ade.md`](./30-sc-pml-ade.md)).

## 2026-09-07 — first pause

Tried Salvador’s volumetric ADE CFS-CPML (Gedney & Zhao), stopped in favor of classical CuDG3D-style \(J\)/\(M\) (later replaced by SC-PML \(P\)-form).

---

## Intent (historical)

| Item | Value |
|------|-------|
| Formulation | Volumetric **CFS-CPML** via **ADE** (auxiliaries \(\psi\)) |
| Reference | Gedney & Zhao, IEEE TAP 2010 |

## Why we stopped

The stretch-derivative driver into \(\psi\) must match Maxwell’s discrete \(D\). Upwind Zero/Two on `globalOperator_` vs centered ADE → closed-loop instability on 2D upwind slabs; the 2026-09-09 centered rebuild still failed 1D absorption.

SC-PML (field-driven \(P\leftarrow F\)) avoids a second discrete \(D\). Do not reintroduce Gedney without a new design review.

---

## Intent (Salvador / Gedney)

| Item | Value |
|------|-------|
| Formulation | Volumetric **CFS-CPML** via **ADE** (auxiliaries \(\psi\) / \(\varphi\)) |
| Reference | Gedney & Zhao, IEEE TAP 2010, DOI [10.1109/TAP.2009.2037765](https://doi.org/10.1109/TAP.2009.2037765) |
| Integration | **`GlobalEvolution` only**; RK4 through `Mult()`; attribute-tagged PML volumes |
| State | Extended ODE `[Ex…Hz \| ψ]`; probes/ParaView export fields only |
| Acceptance | −40 dB reflection (DFT on vacuum probes) |

Region tagging (unchanged for the next formulation): Gmsh volume attributes → JSON `"type": "PML"` + `active_axes` + σ grading. CFS-only JSON fields were `kappa_max` / `alpha_max` (removed from the live schema).

---

## What was built

- Assembled Maxwell curl+flux in `globalOperator_`
- Separate `PMLOperator_`: ADE decay, stretch-derivative driver \(D(F)\), field correction
- 1D **ψ-form** (\(\sigma/\kappa\)-weighted vol+face); dim≥2 **φ-form** (unit SBP \(D\), σ in SPD mass/corr)
- Profiles from mesh depth + JSON grading; multi-axis via `active_axes`
- Diagnostics: operator audit, spectrum probes, hybrid corr experiments (`PML_HYBRID_CORR_THETA`)

Design docs from the CFS era: [`01-physics-and-formulation.md`](./01-physics-and-formulation.md), [`02-codebase-architecture.md`](./02-codebase-architecture.md), [`04-implementation-plan.md`](./04-implementation-plan.md).

---

## How far we got

| Case / check | Result |
|--------------|--------|
| 1D slab DFT (historical `1D_PML_DFT`) | **PASS** (~−40 dB and better; ~−338 dB on the locked DFT case) |
| 2D X-slab **centered** (`upwind_alpha=0`, order ≥ 3) | **PASS** absorption + late-time stable (~−43 to −53 dB) |
| 2D X-slab **upwind** (`upwind_alpha=1`) | Early absorption then **late blow-up** |
| Outer BC (PEC / SMA / `PML_NONE`) | Did **not** drive the 2D spectrum issue |
| Hybrid field corr \(\theta\approx 0.5\) | Stabilized Euclidean spectrum but **killed** absorption (~−9 dB) |
| Marker-based ADE Zero/Two alone | Hurt 1D DFT — discarded |

Primary write-ups: [`21-session-1d-solidification.md`](./21-session-1d-solidification.md), [`23-session-2d-closed-loop-spectrum.md`](./23-session-2d-closed-loop-spectrum.md), [`24-session-2d-hybrid-correction.md`](./24-session-2d-hybrid-correction.md), [`25-session-2d-order-upwind-centered.md`](./25-session-2d-order-upwind-centered.md).

---

## Why we stopped

The hard part of Gedney ADE in DG is that the stretch-derivative driver into \(\psi/\varphi\) must match Maxwell’s discrete \(D\). Global upwind adds Zero/Two on `globalOperator_` that a centered-only ADE driver never saw → closed-loop instability on 2D upwind slabs.

Cudg3d’s classical ADE (\(J\)/\(M\) + volume \(\sigma\)) never builds that stretch-derivative aux driver; PML rides on the **same** curl+flux RHS. Comparison: [`26-cudg3d-pml-comparison.md`](./26-cudg3d-pml-comparison.md).

**Decision (2026-09-07):** remove the CFS ADE / `PMLOperator_` code path; keep region scaffolding (tags, `active_axes`, σ profiles); next implementation is classical volumetric ADE-PML (CuDG3D-style) on MFEM/`GlobalEvolution`.

---

## Code / cases after cleanup

- **Removed:** ψ layout, `buildPMLOperator`, audits/sign tests, CFS coefficients, experimental PML case dirs, `PML_NONE` BC, `kappa_max` / `alpha_max` from JSON.
- **Kept cases:** `1D_PML`, `1D_PML_buffer`, `2D_PML_X_slab`, `2D_RCS_Circle_Vol_PML`.
- **Kept scaffolding:** `PMLProperties`, σ `PMLProfiles`, volume markers, driver `"type": "PML"` parse.
- Until classical ADE is wired, PML-tagged regions behave like vacuum materials (no absorption). Do not use as acceptance runs.
