# Session — 2D X-slab order ↑ and upwind vs centered (2026-09-04)

Follow-up to [`24-session-2d-hybrid-correction.md`](./24-session-2d-hybrid-correction.md).

## Setup

Mesh unchanged (coarse triangles, PML on attr 3, \(x\in[1,2]\)).

Probes:

| Probe | Position | Role |
|-------|----------|------|
| 0 | \((-0.5, 0.5)\) | region 1 — reflected (scattered) |
| 1 | \((0.5, 0.5)\) | region 2 — incident (total-field) |

Code path: dim≥2 **φ-form** ADE, **faces on**, hybrid **off** (θ=0). Global `upwind_alpha` is the only flux switch.

Strange finding recalled from [`21-session-1d-solidification.md`](./21-session-1d-solidification.md): ψ/φ ADE driver is **centered SBP** (vol + OneNormal). Global `upwind_alpha>0` adds Zero/Two on the field curl that the ADE does **not** track → discrete mismatch. Centered (`upwind_alpha=0`) aligns global \(D\) with the ADE driver.

## Results (T=20, Ez, no hybrid)

| Case | order | dt | `upwind_alpha` | \(R_{\mathrm{dB}}\) DFT | late \|Ez\|_ref | max\|ψ\| | Verdict |
|------|-------|-----|----------------|-------------------------|-----------------|---------|---------|
| `2D_PML_X_slab` | 3 | 5e-4 | **1** | −36.8 | grows | \(1.8\times10^5\) @ t≈11 | **FAIL** (late blow) |
| `2D_PML_X_slab_centered` | 3 | 5e-4 | **0** | **−42.8** | \(O(10^{-2})\) | ≈20 | **PASS** absorb + stable |
| upwind | 4 | 2.5e-4 | 1 | −29.1 | grows | \(10^7\) @ t≈9 | **FAIL** |
| centered | 4 | 2.5e-4 | 0 | **−52.5** | \(O(10^{-2})\) | ≈37 | **PASS** (best \(R\)) |

Order-2 φ-form + faces (prior run) already reached early \(R\sim-41\,\mathrm{dB}\) but blew near t≈14–16. Raising order **does not** fix upwind; it **does** lock centered absorption past −40 dB with quiet late fields.

## Takeaway

1. **Use centered slab for 2D PML acceptance** (`2D_PML_X_slab_centered`, order ≥ 3).
2. Upwind global flux remains a late-time hazard on this mesh until the ADE driver shares the same upwind \(D\) as `globalOperator_` (not marker-based Zero/Two).
3. Hybrid \(C\leftrightarrow -B^{\mathsf T}\) still trades absorption for Euclidean stability — not used here.
4. 1D `1D_PML_DFT` still **PASS** (~−338 dB) after these changes.

## Cases left at order 3

- [`2D_PML_X_slab/`](../../testData/maxwellInputs/2D_PML_X_slab/) — upwind stress (known late FAIL)
- [`2D_PML_X_slab_centered/`](../../testData/maxwellInputs/2D_PML_X_slab_centered/) — **preferred** absorb+stable check
