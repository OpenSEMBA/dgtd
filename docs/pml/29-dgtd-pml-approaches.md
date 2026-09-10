# DGTD + PML approaches (survey shortlist)

Literature used to choose the live SC-PML path ([`30-sc-pml-ade.md`](./30-sc-pml-ade.md)).

| Line | Form | Fit for MFEM / GlobalEvolution |
|------|------|--------------------------------|
| Chen/Bagci SC-PML (arXiv:2006.02551) | Field-driven \(P_E,P_H\) + volume tensors \(a,b,c,d(\sigma,\kappa)\) | **Chosen** — same `MassIntegrator` + `AddMult` pattern as classical |
| Lu/Cai/Zhang UPML ADE (JCP 2004) | Polarization currents + Riemann flux with aux | Follow-on if volume SC-PML fails corners |
| Peng et al. UPML DGTD (ISAP 2013) | Aux + Riemann treated dependently | Same as Lu/Cai (face coupling) |
| Gedney & Zhao ADE CFS (TAP 2010) | \(\psi\sim D(F)\) stretch derivative | Parked — discrete-\(D\) mismatch twice ([`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md)) |
| Xiao–Liu well-posed PML | Unsplit strongly hyperbolic | Fallback only |
| Berenger split-field | Split \(E/H\) | Avoid (weakly well-posed / late-time) |

WAA in Bagci is a nodal-DG memory trick; MFEM already assembles QP-graded global masses — take the **PDE**, not WAA.
