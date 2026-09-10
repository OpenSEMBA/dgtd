# Classical ADE-PML (CuDG3D-style) — historical

**Status: replaced** by SC-PML ADE ([`30-sc-pml-ade.md`](./30-sc-pml-ade.md)). At \(\kappa\equiv 1\), Bagci SC-PML with \(J=\sigma^2 P\) is equivalent to the CuDG3D \(J\)/\(M\) system documented here.

## Classical formulation (archive)

Uniaxial or multi-axis stretch (normalized \(\varepsilon=\mu=1\)), volume-only. For each active stretch axis \(s\):

```text
∂t E_s  +=  σ_s E_s − J_s
∂t E_⊥  −=  σ_s E_⊥
∂t J_s   =  σ_s² E_s − σ_s J_s
(same for H / M)
```

Former layout: `[EH | J/M…]` via `ClassicalPMLLayout` (removed). See [`26-cudg3d-pml-comparison.md`](./26-cudg3d-pml-comparison.md).
