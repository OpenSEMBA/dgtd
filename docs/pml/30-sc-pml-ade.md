# SC-PML ADE (Bagci/Chen) — active path

**Status: active** (replaces CuDG3D classical \(J\)/\(M\) and the paused Gedney \(\psi\sim D(F)\) attempts).

## Formulation

Normalized \(\varepsilon=\mu=1\). For component \(u\in\{x,y,z\}\) with cyclic \((u,v,w)\):

```text
∂t (a E)_u = (∇×H)_u − b_uu E_u − c_uu (P_E)_u
∂t (a H)_u = −(∇×E)_u − b_uu H_u − c_uu (P_H)_u
∂t (P_E)_u = κ_u^{-1} E_u − d_uu (P_E)_u
∂t (P_H)_u = κ_u^{-1} H_u − d_uu (P_H)_u

a_uu = κ_v κ_w / κ_u
b_uu = (σ_v κ_w + σ_w κ_v − a_uu σ_u) / κ_u
c_uu = σ_v σ_w − b_uu σ_u
d_uu = σ_u / κ_u
```

At \(\kappa\equiv 1\) (uniaxial), this is equivalent to CuDG3D with \(J=\sigma^2 P\).

### Discrete Mult (MFEM)

1. `globalOperator_` → unit-mass curl+flux into EH.
2. If any `kappa_max > 1`: for each component \(u\), `out_{F_u} += (M_a^{-1}M - M^{-1}M)\, out_{F_u}` (PML-marked), so curl becomes \(M_a^{-1}\) weak form on PML.
3. `scpmlOperator_` AddMult: EH damping with \(M_a^{-1}\mathrm{Mass}(b,c)\) when \(\kappa>1\), else unit \(M^{-1}\); \(P\) rows always unit \(M^{-1}\mathrm{Mass}(1/\kappa,d)\).

## State

```text
[Ex Ey Ez Hx Hy Hz | PEx PEy PEz | PHx PHy PHz]
```

`SCPMLLayout`: \(n_{\mathrm{aux}}=6N\) when any PML region exists.

## JSON

No `pml_formulation` switch. On `"type": "PML"`:

| Field | Default | Notes |
|-------|---------|-------|
| `kappa_max` | `1` | ≥ 1; graded with \(\sigma\); enables curl \(a\)-rescale |
| `alpha_max` | `0` | Must be 0 until CFS pole is wired |
| `active_axes`, `grading_order`, `target_reflection`, `stretch_mode` | as before | Unchanged |

## Acceptance

| Case | Gate |
|------|------|
| `1D_PML` | DFT ≤ −40 dB (\(\kappa=1\)) |
| `1D_PML` + `kappa_max>1` | DFT ≤ −40 dB, finite |
| `2D_PML_X_slab` | Late-stable; centered then upwind |
| Dipole / multi-axis | After 2D upwind |

## References

- Chen et al., arXiv:2006.02551 — SC-PML ADE eqs. 5–9  
- Survey: [`29-dgtd-pml-approaches.md`](./29-dgtd-pml-approaches.md)  
- Gedney archive: [`27-gedney-cfs-paused.md`](./27-gedney-cfs-paused.md)
