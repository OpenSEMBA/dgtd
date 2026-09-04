# 1D_SMA_DFT

SMA-only baseline matched to [`1D_PML_DFT`](../1D_PML_DFT/): same mesh, TFSF, probes, and
solver options; materials are **all vacuum** (no PML). Outer BCs remain SMA at
\(x=-3\) and \(x=3.2\).

Use with [`scripts/pml_dft_reflection.py`](../../../scripts/pml_dft_reflection.py) to
compare reflection against the volumetric PML case.

### 2026-09-04 snapshot (probe 0 at \(x=0\))

| Case | Reflected lobe (time) | \(R_{\mathrm{dB}}(f_{\mathrm{peak}})\) |
|------|-----------------------|----------------------------------------|
| `1D_PML_DFT` | ~0.059 at \(t\approx 9.33\) (PML interface) | **−28.8 dB** (FAIL −40 dB) |
| `1D_SMA_DFT` | numerical noise at SMA return time | **≪ −100 dB** (PASS) |

In 1D normal incidence, first-order SMA is essentially exact for a right-going plane
wave, so this baseline is a very strong absorber. The PML case is stable but still
returns a clear interface reflection above the −40 dB gate — and worse than SMA alone
on this particular 1D test.
