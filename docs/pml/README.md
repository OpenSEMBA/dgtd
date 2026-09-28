# Volumetric PML for OpenSEMBA/dgtd

> **How to mesh and tag a box:** [`designing-a-box-pml.md`](./designing-a-box-pml.md). Absorption reference is the 2D dipole table there.  
> **Cartesian:** SC-PML ADE (Bagci/Chen) — [`30-sc-pml-ade.md`](./30-sc-pml-ade.md).  
> **Radial `"R"`:** removed from the solver; archive [`33-radial-sc-pml-ade.md`](./33-radial-sc-pml-ade.md). Parser rejects `active_axes: ["R"]`.  
> **Status:** [`31-sc-pml-status.md`](./31-sc-pml-status.md). Gedney CFS paused — [`27`](./27-gedney-cfs-paused.md).

## Document index

| File | Purpose |
|------|---------|
| [designing-a-box-pml.md](./designing-a-box-pml.md) | Mesh cut, JSON attributes, dipole absorption table |
| [33-radial-sc-pml-ade.md](./33-radial-sc-pml-ade.md) | Archive: former `"R"` Cyl UPML |
| [32-radial-dgtd-literature.md](./32-radial-dgtd-literature.md) | Radial literature notes |
| [31-sc-pml-status.md](./31-sc-pml-status.md) | Status report |
| [30-sc-pml-ade.md](./30-sc-pml-ade.md) | Cartesian SC-PML ADE |
| [29-dgtd-pml-approaches.md](./29-dgtd-pml-approaches.md) | Literature shortlist |
| [27-gedney-cfs-paused.md](./27-gedney-cfs-paused.md) | Gedney CFS archive |

**Cases:** Cartesian slabs and the 2D dipole (`2D_Dipole_PML`, `2D_Dipole_PML_GO*`). Radial PML cases are removed.

## Solver path

Only **GlobalEvolution** is in scope.
