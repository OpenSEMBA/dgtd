# Locked decisions (do not re-litigate without user approval)

**Updated 2026-09-09:** Active PML is **SC-PML ADE** (Bagci/Chen). Gedney \(\psi\sim D(F)\) remains archived. Classical CuDG3D \(J/M\) replaced by equivalent \(\kappa\equiv 1\) SC-PML \(P\)-form.

## Formulation

| Decision | Value |
|----------|-------|
| PML type (**active**) | Volumetric **SC-PML ADE** (\(P_E/P_H\) + volume \(a,b,c,d(\sigma,\kappa)\)) |
| Primary reference | Chen et al., arXiv:2006.02551 (SC-PML ADE) |
| Gedney CFS | **Parked** — [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md) |
| Surface PML | **Rejected** — no `SBC_PML` |

## Solver integration

| Decision | Value |
|----------|-------|
| Evolution operator | **`GlobalEvolution` only** |
| Time integrator | **RK4** (explicit), PML through **`Mult()`** |
| `ImplicitSolve()` | **Deferred** |
| Coexistence | **SGBC** and **TFSF** remain |

## Mesh and JSON

| Decision | Value |
|----------|-------|
| PML region definition | **Gmsh volume → attribute tag → JSON `"type": "PML"`** |
| Thickness | **From mesh geometry** |
| Vacuum match | **`matches_vacuum: true`** |
| `bulk_conductivity` on PML | **Forbidden** |
| Grading | `grading_order`, `target_reflection`, `active_axes`, `kappa_max` (≥1) |
| `pml_formulation` switch | **Removed** |
| `alpha_max` | Must be **0** until CFS pole is wired |

## State vector and I/O

| Decision | Value |
|----------|-------|
| ODE state | `[E/H (6N); P_E (3N); P_H (3N)]` via `SCPMLLayout` |
| Probes / Paraview / MOR | **E and H only** |
| Units | **Normalized** (`c = ε = μ = 1`) |

## Acceptance

| Decision | Value |
|----------|-------|
| Reflection target | **−40 dB** DFT on vacuum probes |
| Kept reference cases | `1D_PML`, `1D_PML_buffer`, `2D_PML_X_slab`, `2D_RCS_Circle_Vol_PML` |

## What remains rejected

1. Spatially varying **`bulk_conductivity`** as a fake PML.
2. SGBC-style **sub-solver** for volumetric PML.
3. Reintroducing **Gedney \(\psi\sim D(F)\)** without a new discrete-\(D\) design review.
4. A **1D-only** PML module fork.
