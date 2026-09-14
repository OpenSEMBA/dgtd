# Volumetric PML for OpenSEMBA/dgtd

> **Cartesian:** SC-PML ADE (Bagci/Chen) — [`30-sc-pml-ade.md`](./30-sc-pml-ade.md).  
> **Radial `"R"`:** removed from the solver; archive [`33-radial-sc-pml-ade.md`](./33-radial-sc-pml-ade.md). Parser rejects `active_axes: ["R"]`.  
> **Status:** [`31-sc-pml-status.md`](./31-sc-pml-status.md). Gedney CFS paused — [`27`](./27-gedney-cfs-paused.md).

## Document index

| File | Purpose |
|------|---------|
| [33-radial-sc-pml-ade.md](./33-radial-sc-pml-ade.md) | Archive: former `"R"` Cyl UPML |
| [32-radial-dgtd-literature.md](./32-radial-dgtd-literature.md) | Radial literature notes |
| [31-sc-pml-status.md](./31-sc-pml-status.md) | Status report |
| [30-sc-pml-ade.md](./30-sc-pml-ade.md) | Cartesian SC-PML ADE |
| [29-dgtd-pml-approaches.md](./29-dgtd-pml-approaches.md) | Literature shortlist |
| [27-gedney-cfs-paused.md](./27-gedney-cfs-paused.md) | Gedney CFS archive |

**Cases:** Cartesian slabs/dipole. Onion `"R"` JSON under `2D_Onion_*_R` / `*PML_R*` stays on disk as expected parse failures.

## Solver path

Only **GlobalEvolution** is in scope.
