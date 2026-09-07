# Locked decisions (do not re-litigate without user approval)

This document records decisions from the planning conversation. **Updated 2026-09-07:** Gedney ADE CFS-CPML is **paused**; classical volumetric ADE-PML (CuDG3D-style) is the next active path. Archive: [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md).

## Formulation

| Decision | Value |
|----------|-------|
| PML type (paused) | Volumetric **CFS-CPML** via ADE \(\psi\) — **not active in code** |
| PML type (**next / active**) | Volumetric **classical ADE-PML** (\(J\)/\(M\) + volume \(\sigma\)), CuDG3D-style on MFEM |
| CFS primary reference | Gedney & Zhao, IEEE TAP 2010, DOI `10.1109/TAP.2009.2037765` (historical) |
| Secondary reference | Taflove for stretch / \(\sigma\) grading intuition — **not** FDTD stencils |
| Surface PML | **Rejected** — no `SBC_PML` |

## Solver integration

| Decision | Value |
|----------|-------|
| Evolution operator | **`GlobalEvolution` only** |
| Time integrator | **RK4** (explicit), PML dynamics through **`Mult()`** |
| `ImplicitSolve()` | **Deferred** until explicit validation |
| Regional IMEX | **Out of scope** for v1 |
| Coexistence | **SGBC** and **TFSF** remain; not applied on PML surfaces in user workflows |

## Dimensionality

| Decision | Value |
|----------|-------|
| Code structure | **Dimension-agnostic** (loop X/Y/Z, skip `d >= mesh.Dimension()`) |
| Multi-axis | Per-tag `active_axes` (uniaxial slabs + multi-block RCS) |

## Mesh and JSON (retained)

| Decision | Value |
|----------|-------|
| PML region definition | **Gmsh volume → attribute tag → JSON `"type": "PML"`** |
| Thickness | **From mesh geometry** — no `pml_thickness` JSON |
| Vacuum match | **`matches_vacuum: true`** → ε = μ = 1 |
| `bulk_conductivity` on PML | **Forbidden** |
| Grading | `grading_order`, `target_reflection`, `active_axes` |
| CFS-only fields | **`kappa_max` / `alpha_max` removed** from schema and parser |

## State vector and I/O

| Decision | Value |
|----------|-------|
| ODE state (CFS, paused) | Was `[E; ψ]` — removed from code |
| ODE state (**classical, active**) | `[E/H (6N); J/M aux]` via `ClassicalPMLLayout` |
| Probes / Paraview / MOR | **E and H only**; auxiliaries not exported |
| Units | **Normalized** (`c = ε = μ = 1`) |

## Acceptance

| Decision | Value |
|----------|-------|
| Probes | In **vacuum** |
| Frequency analysis | **DFT offline** |
| Reflection target | **−40 dB** |
| Kept reference cases | `1D_PML`, `1D_PML_buffer`, `2D_PML_X_slab`, `2D_RCS_Circle_Vol_PML` |

## What remains rejected

1. Spatially varying **`bulk_conductivity`** as a fake PML.
2. SGBC-style **sub-solver** for volumetric PML.
3. FDTD-only / Yee-tied path.
4. Keeping **`SBC_PML`** alongside volumetric PML.
5. A **1D-only** PML module fork.
