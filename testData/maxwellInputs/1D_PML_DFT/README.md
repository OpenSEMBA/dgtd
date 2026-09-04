# 1D_PML_DFT

Vacuum-probe DFT acceptance case for volumetric CFS-CPML.

- **Mesh:** copy of `1D_PML_buffer` (PML on \(x\in[2,3]\), vacuum buffer to \(3.2\), TFSF at \(x=-2.5\)).
- **Probes (vacuum only):**

| ID | x | Role |
|----|---|------|
| 0 | 0.0 | Primary DFT (best lobe separation) |
| 1 | 1.0 | Check |
| 2 | 1.5 | Check (closer to interface) |

### Timing (code time, `dt=0.001`)

| Probe | Incident peak | Reflected peak |
|-------|---------------|----------------|
| x=0 | ~5.33 | ~9.33 |
| x=1 | ~6.33 | ~8.33 |
| x=1.5 | ~6.83 | ~8.05 |

Recommended DFT windows for probe 0:

```text
--inc-window 3.0 7.2 --ref-window 8.5 11.5
```

See [scripts/README.md](../../../scripts/README.md) and run with `scripts/pml_dft_reflection.py`.
