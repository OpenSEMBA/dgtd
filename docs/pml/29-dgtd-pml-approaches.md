# DGTD + PML approaches (survey shortlist)

| Line | Form | Fit for MFEM / GlobalEvolution |
|------|------|--------------------------------|
| Chen/Bagci SC-PML (arXiv:2006.02551) | Field-driven $P_E,P_H$ + $a,b,c,d(\sigma,\kappa)$ | **Live for Cartesian** `X`/`Y`/`Z` — [`30-sc-pml-ade.md`](./30-sc-pml-ade.md) |
| Ji/Lu/Cai UPML ADE + face aux + Teixeira | Polarization $P,Q$ + cylindrical metric + aux in flux | **Removed** — archive [`33-radial-sc-pml-ade.md`](./33-radial-sc-pml-ade.md) |
| Bagci SC-PML + Teixeira (rotated) | Same ADE in $(r,\theta)$ + $R\,T'\,R^{\mathsf{T}}$ | **Removed** (remanents / bounce) |
| Ji/Lu/Cai volume-only graft | UPML ADE, vacuum Riemann | **Removed** (growth) |
| Peng et al. UPML DGTD (ISAP 2013) | Aux + Riemann treated dependently | Face follow-on (Phase 4b) |
| Kapidani/Schöberl stretch maps | Covariant complex stretch | Parked — matrix-free / map-based |
| Duru et al. energy-stable | Aux in face fluxes | Phase 4b fallback |
| Gedney & Zhao ADE CFS (TAP 2010) | $\psi\sim D(F)$ stretch derivative | Parked — [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md) |
| Berenger split-field | Split $E/H$ | Avoid |

WAA in Bagci is a nodal-DG memory trick; MFEM already assembles QP-graded global masses — take the **PDE**, not WAA.
