# Volumetric PML for OpenSEMBA/dgtd

> **Active path:** classical volumetric ADE-PML (CuDG3D-style \(J\)/\(M\)) — [`28-classical-ade-pml.md`](./28-classical-ade-pml.md).  
> **CFS-CPML is paused** (2026-09-07) — [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md). Session docs 09–26 are historical.

Static design memory for PML work in **GlobalEvolution**. Region tagging and σ grading remain live scaffolding; ADE absorption operators are not.

## Document index

| File | Purpose |
|------|---------|
| [00-decisions-locked.md](./00-decisions-locked.md) | Product decisions (updated for CFS pause → classical ADE next) |
| [01-physics-and-formulation.md](./01-physics-and-formulation.md) | *(historical)* CFS-CPML physics, stretch maps, ADE structure |
| [02-codebase-architecture.md](./02-codebase-architecture.md) | *(historical)* How CFS PML fit into `GlobalEvolution` |
| [03-json-schema-and-mesh.md](./03-json-schema-and-mesh.md) | JSON/Gmsh contract (σ grading; no κ/α) |
| [04-implementation-plan.md](./04-implementation-plan.md) | *(historical)* CFS milestone plan |
| [05-verification.md](./05-verification.md) | Acceptance tests, −40 dB criterion, probe workflow |
| [06-cleanup-inventory.md](./06-cleanup-inventory.md) | *(historical)* SBC_PML / leftover cleanup |
| [07-time-integration.md](./07-time-integration.md) | RK4 / `Mult()` vs `ImplicitSolve()`, IMEX deferral |
| [08-conversation-snapshot.md](./08-conversation-snapshot.md) | Planning decisions from initial design conversation |
| [09-session-record-2026-05-27.md](./09-session-record-2026-05-27.md) | **2026-05-27** CFS implementation log |
| [10-session-2026-05-28-stability-dft.md](./10-session-2026-05-28-stability-dft.md) | **2026-05-28** t=20 stability, α sweep |
| [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md) | **2026-05-28** Sign sweep S1–S7 |
| [12-session-pml-boundary-terminator.md](./12-session-pml-boundary-terminator.md) | **2026-05-28** `PML_NONE` experiments |
| [13-session-centered-flux-and-deriv-split.md](./13-session-centered-flux-and-deriv-split.md) | **2026-05-28** centered flux split (reverted) |
| [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md) | **2026-05-27** Revert + audit |
| [15-session-operator-duplicate-audit.md](./15-session-operator-duplicate-audit.md) | **2026-05-27** Operator duplicate audit |
| [17-session-centered-frame-probe-correlation.md](./17-session-centered-frame-probe-correlation.md) | **2026-05-27** Frame/probe correlation |
| [18-session-centered-vs-upwind-correlation.md](./18-session-centered-vs-upwind-correlation.md) | **2026-05-27** Centered vs upwind |
| [19-session-sign-matrix-interface-window.md](./19-session-sign-matrix-interface-window.md) | **2026-05-27** Sign matrix at t=8.05 |
| [20-pipeline-reexploration.md](./20-pipeline-reexploration.md) | **2026-09-04** Full pipeline audit |
| [21-session-1d-solidification.md](./21-session-1d-solidification.md) | **2026-09-04** 1D CFS findings (PASS) |
| [22-session-2d-assembly-instability.md](./22-session-2d-assembly-instability.md) | **2026-09-04** 2D assembly instability |
| [23-session-2d-closed-loop-spectrum.md](./23-session-2d-closed-loop-spectrum.md) | **2026-09-04** driver↔corr spectrum |
| [24-session-2d-hybrid-correction.md](./24-session-2d-hybrid-correction.md) | **2026-09-04** Hybrid corr experiment |
| [25-session-2d-order-upwind-centered.md](./25-session-2d-order-upwind-centered.md) | **2026-09-04** Centered PASS; upwind FAIL |
| [26-cudg3d-pml-comparison.md](./26-cudg3d-pml-comparison.md) | **2026-09-07** Cudg3d J/M vs CFS ψ |
| [27-gedney-cfs-paused.md](./27-gedney-cfs-paused.md) | **2026-09-07** CFS pause archive |
| [28-classical-ade-pml.md](./28-classical-ade-pml.md) | **2026-09-07** Classical ADE J/M wiring (active) |

**DFT tooling:** [`scripts/pml_dft_reflection.py`](../../scripts/pml_dft_reflection.py) (still useful for classical ADE acceptance).  
**Kept cases:** [`1D_PML`](../../testData/maxwellInputs/1D_PML/), [`1D_PML_buffer`](../../testData/maxwellInputs/1D_PML_buffer/), [`2D_PML_X_slab`](../../testData/maxwellInputs/2D_PML_X_slab/), [`2D_RCS_Circle_Vol_PML`](../../testData/maxwellInputs/2D_RCS_Circle_Vol_PML/).

## Primary references

1. **Gedney, S. D. and Zhao, B.**, *An Auxiliary Differential Equation Formulation for the Complex-Frequency Shifted PML*, IEEE Trans. Antennas Propag., vol. 58, no. 3, pp. 838–847, Mar. 2010.  
   **DOI:** [10.1109/TAP.2009.2037765](https://doi.org/10.1109/TAP.2009.2037765) — *CFS attempt reference; paused.*

2. **Classical ADE-PML (next):** OpenSEMBA Cudg3d polarization-current ADE (\(J\)/\(M\) + volume \(\sigma\)) — see [`26-cudg3d-pml-comparison.md`](./26-cudg3d-pml-comparison.md) and [`external/Cudg3d/`](../../external/Cudg3d/).

3. **Institutional context:** Salvador González García (UGR / OpenSEMBA / GEG). Surface `SBC_PML` remains rejected.

## One-line goal (current)

Attribute-tagged **volumetric classical ADE-PML** in **GlobalEvolution** (wired — see [`28-classical-ade-pml.md`](./28-classical-ade-pml.md)), validated on kept 1D/2D cases (−40 dB DFT).

## Solver path

Only **GlobalEvolution** is in scope. **MaxwellEvolution** and **HesthavenEvolution** are deprecated for new features.
