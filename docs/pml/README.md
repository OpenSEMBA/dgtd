# Volumetric PML for OpenSEMBA/dgtd

> **Active path:** SC-PML ADE (Bagci/Chen \(P_E/P_H\)) — [`30-sc-pml-ade.md`](./30-sc-pml-ade.md).  
> **Status report (math, impl, dipole GO vs SMA):** [`31-sc-pml-status.md`](./31-sc-pml-status.md).  
> **Gedney CFS paused** — [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md). Approaches: [`29-dgtd-pml-approaches.md`](./29-dgtd-pml-approaches.md).

Static design memory for PML work in **GlobalEvolution**.

## Document index

| File | Purpose |
|------|---------|
| [00-decisions-locked.md](./00-decisions-locked.md) | Product decisions |
| [31-sc-pml-status.md](./31-sc-pml-status.md) | **Status report** — physics, ADE vs Gedney, Mult wiring, GO vs SMA |
| [30-sc-pml-ade.md](./30-sc-pml-ade.md) | Live SC-PML ADE equations / JSON |
| [29-dgtd-pml-approaches.md](./29-dgtd-pml-approaches.md) | Literature shortlist |
| [27-gedney-cfs-paused.md](./27-gedney-cfs-paused.md) | Gedney CFS archive |
| [28-classical-ade-pml.md](./28-classical-ade-pml.md) | Historical CuDG3D \(J/M\) (replaced by 30) |

**DFT tooling:** [`scripts/pml_dft_reflection.py`](../../scripts/pml_dft_reflection.py), [`scripts/dipole_pml_grading_dft.py`](../../scripts/dipole_pml_grading_dft.py).  
**Kept / sweep cases:** [`1D_PML`](../../testData/maxwellInputs/1D_PML/), [`2D_PML_X_slab`](../../testData/maxwellInputs/2D_PML_X_slab/), [`2D_Dipole`](../../testData/maxwellInputs/2D_Dipole/), [`2D_Dipole_PML_GO*`](../../testData/maxwellInputs/2D_Dipole_PML_GO0/).

## Primary references

1. **Chen et al.**, arXiv:2006.02551 — live SC-PML ADE.  
2. **Gedney & Zhao**, IEEE TAP 2010 — CFS \(\psi\) form; *paused* (DOI [10.1109/TAP.2009.2037765](https://doi.org/10.1109/TAP.2009.2037765)).  
3. OpenSEMBA Cudg3d \(J/M\) — \(\kappa\equiv 1\) limit; [`26-cudg3d-pml-comparison.md`](./26-cudg3d-pml-comparison.md).

## One-line goal (current)

Attribute-tagged **volumetric SC-PML ADE** in **GlobalEvolution**, validated on slabs (−40 dB) and dipole GO vs SMA (~−70 dB late reflection vs ~−20 dB SMA).

## Solver path

Only **GlobalEvolution** is in scope. **MaxwellEvolution** and **HesthavenEvolution** are deprecated for new features.
