# Volumetric CFS-CPML for OpenSEMBA/dgtd

Static design memory for implementing **Complex-Frequency-Shifted Convolutional PML (CFS-CPML)** as **Auxiliary Differential Equations (ADE)** in the **GlobalEvolution** DGTD solver.

This folder exists so future Agent sessions and human reviewers do not depend on chat context. Read documents in numeric order for a full picture.

## Document index

| File | Purpose |
|------|---------|
| [00-decisions-locked.md](./00-decisions-locked.md) | Non-negotiable product and scope decisions agreed with the user |
| [01-physics-and-formulation.md](./01-physics-and-formulation.md) | CFS-CPML physics, references, stretch maps, ADE structure, Maxwell mapping |
| [02-codebase-architecture.md](./02-codebase-architecture.md) | How PML fits into `GlobalEvolution`, `Fields`, `DGOperatorFactory`, etc. |
| [03-json-schema-and-mesh.md](./03-json-schema-and-mesh.md) | JSON contract, Gmsh tagging, depth profiles, no separate thickness JSON |
| [04-implementation-plan.md](./04-implementation-plan.md) | Milestones, file touch list, implementation order, branch workflow |
| [05-verification.md](./05-verification.md) | Acceptance tests, −40 dB criterion, probe workflow |
| [06-cleanup-inventory.md](./06-cleanup-inventory.md) | Code to remove (`SBC_PML`, leftovers), code to keep |
| [07-time-integration.md](./07-time-integration.md) | RK4 / `Mult()` vs `ImplicitSolve()`, IMEX deferral |
| [08-conversation-snapshot.md](./08-conversation-snapshot.md) | Planning decisions from initial design conversation |
| [09-session-record-2026-05-27.md](./09-session-record-2026-05-27.md) | **2026-05-27** implementation log: factory `PMLOperator_`, face driver fix, vector ψ |
| [10-session-2026-05-28-stability-dft.md](./10-session-2026-05-28-stability-dft.md) | **2026-05-28** t=20 stability, α sweep, probe/DFT windows |
| [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md) | **2026-05-28** Sign sweep S1–S7, outer-PML diagnostics, split-audit notes |
| [12-session-pml-boundary-terminator.md](./12-session-pml-boundary-terminator.md) | **2026-05-28** `PML_NONE` + vacuum-buffer outer BC experiments |
| [13-session-centered-flux-and-deriv-split.md](./13-session-centered-flux-and-deriv-split.md) | **2026-05-28** `upwind_alpha`, centered cases (split experiment — reverted) |
| [14-session-pml-revert-and-audit.md](./14-session-pml-revert-and-audit.md) | **2026-05-27** Revert split, PML ψ upwind driver, probe gate, t=20 matrix |
| [15-session-operator-duplicate-audit.md](./15-session-operator-duplicate-audit.md) | **2026-05-27** Operator assembly / duplicate-curl audit (`PML_OPERATOR_AUDIT`) |
| [17-session-centered-frame-probe-correlation.md](./17-session-centered-frame-probe-correlation.md) | **2026-05-27** Centered case: ParaView frames 48/61/68–75 vs probes + audit |
| [18-session-centered-vs-upwind-correlation.md](./18-session-centered-vs-upwind-correlation.md) | **2026-05-27** ψ Mult diag + centered vs upwind at frames 48–75 |
| [19-session-sign-matrix-interface-window.md](./19-session-sign-matrix-interface-window.md) | **2026-05-27** S0–S6 at t=8.05 on `1D_PML_centered` |

## Primary references

1. **Gedney, S. D. and Zhao, B.**, *An Auxiliary Differential Equation Formulation for the Complex-Frequency Shifted PML*, IEEE Trans. Antennas Propag., vol. 58, no. 3, pp. 838–847, Mar. 2010.  
   **DOI:** [10.1109/TAP.2009.2037765](https://doi.org/10.1109/TAP.2009.2037765)

2. **Secondary sanity check:** Taflove, *Computational Electrodynamics: The Finite-Difference Time-Domain Method* — CPML stretch profiles and implementation intuition (adapt to DGTD, do not copy FDTD stencils).

3. **Institutional context:** Salvador González García (UGR / OpenSEMBA / GEG). Formulation choice is ADE CFS-PML, not FDTD split-field Berenger, not surface `SBC_PML`.

## One-line goal

Add **attribute-tagged volumetric CFS-CPML** to **GlobalEvolution** via **extended ODE state** `[E, H, ψ]`, **dimension-agnostic** from day one, validated first on a **1D mesh** the user will provide, structurally aligned with **`testData/maxwellInputs/2D_RCS_Circle_Vol_PML/`**.

## Explicit non-goals (Milestone A)

- GPU (`GlobalEvolution.cu`) PML kernels
- MPI-specific PML communication beyond existing field exchange
- Regional IMEX (explicit vacuum / implicit PML only)
- Automated gtest cases in `test/cases/` (manual `testData/maxwellInputs` only)
- Paraview export of auxiliary fields ψ
- Surface absorber shortcut `SBC_PML`
- `VolumetricRegionSubMesher` as a runtime dependency (optional diagnostic only)

## Solver path

Only **GlobalEvolution** is in scope. **MaxwellEvolution** and **HesthavenEvolution** are deprecated for new features.
