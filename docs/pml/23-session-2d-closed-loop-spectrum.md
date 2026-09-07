# Session — 2D closed-loop spectrum (2026-09-04)

Follow-up to [`22-session-2d-assembly-instability.md`](./22-session-2d-assembly-instability.md).

## Question

Is the 2D X-slab blow-up from outer BCs, or from the PML operator spectrum?

## Answers (CSR eigenspectra)

### Wall BCs

| Operator | PEC / SMA / `PML_NONE` |
|----------|-------------------------|
| `pml_full.csr` | **Bit-identical** |
| `global_full.csr` | max Re ≈ 0 (PEC/SMA) or 0.04 (`PML_NONE`); max \|λ\| ≈ 61 |

Outer BC swap changes late growth rate but **not** the PML spectrum. Walls are not the eigenvalue source.

### Block ablation (`pml_full`, before further experiments)

| Blocks | Re max |
|--------|--------|
| mass / driver / corr alone | **≤ 0** |
| mass+driver, mass+corr | **≤ 0** |
| **driver+corr** | **+16** |
| **full** | **+7.3** |

Instability is the **driver ↔ field-correction closed loop**. Unstable eigenvectors are ψ-dominated.

### 1D reference

`1D_PML_buffer` `pml_full`: Re max **≤ 0**.

## Attempts that did **not** stabilize Re(λ)

1. Nodal σ row-scale of unweighted / global-style D (matches nodal decay collocated form).
2. Lumped `M⁻¹` + nodal σ (Euclidean product for D).
3. Weak Galerkin corr with nodal decay.
4. `PMLElementConstCoefficient` + full Galerkin decay/driver/corr — mass-only still Re>0 because **global `M⁻¹` × PML-only `Mass`** is not `c I`.
5. Euclidean skew-symmetrize of vol/face driver blocks.
6. Stronger nodal decay (even ×10 leaves Re>0).

Spectral probe: replacing corr `C` with **`-Bᵀ`** (discrete adjoint of driver) yields Re≈0, but that changes the PDE (not Gedney `u̇ ±= ψ`). It shows the current `(B, C=±I)` pair is not a dissipative structure for this discrete `B`.

## Root-cause picture

- Open-loop nodal decay is fine (Re≤0).
- Discrete driver `B ≈ σ D` on this mesh has a **large Euclidean-symmetric part** (variable σ + DG faces); with nodal `C ≈ ±I` the block `[0, C; B, −Σ]` has positive-Re eigenvalues.
- 1D works because decay, driver, and corr are a **consistent Galerkin** triple (`M⁻¹ Mass` / `M⁻¹(σD)`).

## Code left in tree

Superseded by [`24-session-2d-hybrid-correction.md`](./24-session-2d-hybrid-correction.md): element-const Galerkin + volume-only φ-form + hybrid field correction (θ=0.5) clears Re≤0 on this slab.

## Next (at time of this note; see session 24)

1. Shared discrete `D` with global curl **in the same inner product** as decay/corr (full Galerkin on a subspace, or true collocation D).
2. Or redesign corr as the **M-adjoint** of the driver while preserving Gedney `1/s` meaning (not raw `-Bᵀ` in Euclidean coordinates).
3. Finer PML mesh / milder intra-element σ as a **mitigation**, not a substitute for (1)/(2).
