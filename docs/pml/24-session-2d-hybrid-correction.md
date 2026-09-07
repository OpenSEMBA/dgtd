# Session — 2D Galerkin ADE + hybrid field correction (2026-09-04)

Follow-up to [`23-session-2d-closed-loop-spectrum.md`](./23-session-2d-closed-loop-spectrum.md).

## Goal

Stabilize the 2D ADE closed loop (`driver ↔ field correction`) so `pml_full` has Re(λ)≤0 and `2D_PML_X_slab` runs to T=20 without ψ blow-up, without breaking 1D.

## What landed (dim ≥ 2)

1. **Element-const Galerkin decay / σ / 1/κ** (`PMLElementConstCoefficient`) with matched mass IR (`2p+3`) so `M⁻¹ Mass(c)` stays SPD (mass-only Re≤0).
2. **φ-form ADE**: unit stretch derivative `D` drives φ; σ lives in SPD mass terms (decay + field corr `Mass(σ/κ)`).
3. **Volume-only driver by default** (SBP faces → Re≈6–8). `PML_FORCE_FACE_DRIVER=1` restores faces.
4. **Hybrid field correction** (default θ=0.5):
   \[
   C_h = (1-\theta)\,C_{\mathrm{Gedney}} + \theta\,(-B^{\mathsf T})
   \]
   Applied after CSR merge when `nbrDofs==0`. Override / disable: `PML_HYBRID_CORR_THETA` (`0` = pure Gedney).

### 1D unchanged path

- ψ-form with σ-weighted vol+face, Maxwell `M⁻¹` for field corr, default MFEM IR for scalar inverse mass, face σ = `max(σ+,σ−)`.

## Spectra (`2D_PML_X_slab`, volume-only)

| Config | Re max |
|--------|--------|
| Pure Gedney (θ=0) | ≈ +0.11 |
| Hybrid θ=0.5 (default) | ≈ 1e−5 |
| Pure adjoint θ=1 | ≈ 0 |

## Time domain

| Case | Result |
|------|--------|
| `2D_PML_X_slab` T=8 | peak \|ψ\| ≈ 0.83, final ≈ 0.61 |
| `2D_PML_X_slab` T=20 | peak \|ψ\| ≈ 1.16, late ≥15 bounded, final ≈ 0.52 |
| `1D_PML_buffer` T=20 | peak \|ψ\| ≈ 28, final ~1e−13 (matches prior healthy run) |

## Honest caveat

θ>0 is **not** pure Gedney `u̇ ∓= ψ/κ`. It mixes in the Euclidean adjoint of the discrete driver so the `(B,C)` block cannot inject energy. Empirically θ=0.5 keeps absorption (ψ stays O(1) and decays after the pulse) while clearing Re>0. A future pure-Gedney fix still wants a discrete `D` that is skew in the same product as `C≈±I` (true SBP with consistent σ), so θ can return to 0.

## Env knobs

| Env | Role |
|-----|------|
| `PML_HYBRID_CORR_THETA` | Hybrid weight (default 0.5 on dim≥2; `0` disables) |
| `PML_FORCE_FACE_DRIVER=1` | Include SBP faces in 2D driver |
| `PML_SKIP_FACE_DRIVER=1` | Force volume-only (also 1D) |
| `PML_DECAY_OVERDRIVE` | Scale ADE decay (experiment) |
| `PML_SKEW_DRIVER=1` | Skew-symmetrize weak `D` (worsened spectrum here) |

## Code touch points

- `src/components/DGOperatorFactory.h` — assembly, hybrid post-process
- `src/components/PMLCoefficients.*` — element-const profiles
- Face σ remains `max` (arithmetic mean broke 1D)
