# Locked decisions (do not re-litigate without user approval)

This document records decisions from the planning conversation. Implementation must follow these unless the user explicitly changes them.

## Formulation

| Decision | Value |
|----------|-------|
| PML type | **Volumetric CFS-CPML** via **ADE** (auxiliary differential equations) |
| Primary reference | Gedney & Zhao, IEEE TAP 2010, DOI `10.1109/TAP.2009.2037765` |
| Secondary reference | Taflove CPML (book) for stretch intuition only — **not** FDTD discretization |
| Naming in docs/JSON | **CFS-CPML** / **ADE CFS-PML** (Salvador's vocabulary) |
| UPML | Same physics family; use Gedney ADE paper as source of truth, not a competing implementation |

## Solver integration

| Decision | Value |
|----------|-------|
| Evolution operator | **`GlobalEvolution` only** |
| Time integrator (Milestone A) | **RK4** (explicit), all dynamics through **`Mult()`** including ψ |
| `ImplicitSolve()` | **Deferred** until explicit validation; same operator split, not regional IMEX |
| Regional IMEX (explicit vacuum / implicit PML) | **Out of scope** for v1 |
| Coexistence | **SGBC** and **TFSF** remain; they are never applied on PML surfaces/interfaces in user workflows |

## Dimensionality

| Decision | Value |
|----------|-------|
| Code structure | **Dimension-agnostic** from day one (same pattern as existing `DGOperatorFactory`: loop X/Y/Z, skip `d >= mesh.Dimension()`) |
| 1D test mesh | **Validation only** — not a separate code path |
| Multi-axis corners | Handled by generic formulation + per-tag `active_axes` — not a later "Phase 4 code fork" |

## Mesh and JSON

| Decision | Value |
|----------|-------|
| PML region definition | **Gmsh volume/surface → element attribute tag → JSON `"type": "PML"`** |
| Physical layer thickness | **From mesh geometry** (extent of tagged elements) — **no required `pml_thickness` JSON** |
| `VolumetricRegionSubMesher` | **Not required** for Milestone A; attribute + neighbor logic suffices |
| Vacuum match | **`matches_vacuum: true`** → ε = μ = 1 in PML; no bulk Ohmic loss in PML tags |
| `bulk_conductivity` in PML | **Forbidden / ignored** — PML σ is stretch conductivity, not `Material::sigma_` |

## State vector and I/O

| Decision | Value |
|----------|-------|
| ODE state | `[E (6×N); ψ (n_aux)]` when PML present; `n_aux = 0` when no PML tags |
| Probes / Paraview / MOR | **E and H only** (first `6×ndofs` block); ψ never exported |
| Units | **Normalized** (`c = ε = μ = 1`, same as rest of dgtd) |

## Cleanup and git

| Decision | Value |
|----------|-------|
| Surface PML | **Remove `SBC_PML`** and related driver/README docs |
| Old vol PML experiments | **Remove dead code**; leftover JSON fields were placeholders |
| Branch | **Delete and recreate `vol_pml` from current `dev`** before implementation |
| Prior branch content | **Do not merge** old `vol_pml` attempt — start fresh |

## Platform scope (Milestone A)

| In scope | Out of scope |
|----------|--------------|
| CPU, serial, single rank | GPU PML in `GlobalEvolution.cu` |
| Generic 1D/2D/3D PML operators | Dimension-specific evolvers |
| Manual cases under `testData/maxwellInputs/` | New automated gtests in `test/cases/` |

## Boundaries

| Decision | Value |
|----------|-------|
| Outer domain face terminating PML | Use **SMA** (or existing best absorber) on outer mesh boundary — user test PEC outers were exploratory |
| Vacuum–PML interface | **Standard DG flux**; stretch profile zero at interface (κ→1, σ→0, α→0) |
| PEC on real conductors | Unchanged; separate from PML outer termination policy |

## Acceptance

| Decision | Value |
|----------|-------|
| Probes | User places probes in **vacuum** region |
| Frequency analysis | User performs **DFT offline** on probe time series |
| Reflection target | **−40 dB** (reflected vs incident, frequency domain) |
| Reference case structure | `testData/maxwellInputs/2D_RCS_Circle_Vol_PML/2D_RCS_Circle_Vol_PML.json` |
| First validation mesh | User-provided **1D slab** with PML on one side (e.g. right) |

## What was explicitly rejected

1. Implementing PML as spatially varying **`bulk_conductivity`** in `buildSigmaPiecewiseVector()`.
2. SGBC-style **sub-solver** for volumetric PML (wrong RK stage semantics for ψ).
3. FDTD-only implementation path tied to Yee stencils.
4. Keeping **`SBC_PML`** as parallel feature to volumetric PML.
5. "1D-only PML module" as a shortcut — violates generalist architecture.
