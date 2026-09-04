# 1D_PML_L2

Thickened volumetric PML acceptance case. Same layout as `1D_PML_DFT` / buffer, but
PML thickness **L = 2** on \(x\in[2,4]\) (20 elements at \(\Delta x=0.1\)), vacuum buffer
to \(x=4.2\).

Tags (unchanged semantics):

| Tag | Region |
|-----|--------|
| 1 | Vacuum ends (left \([-3,-2.5]\), right buffer past PML) |
| 2 | Main vacuum \([-2.5, 2]\) |
| 3 | PML |
| Point 1 / 2 | SMA left / right |
| Point 3 | TFSF |
| Point 4 | Vacuum–PML interface at \(x=2\) |

Default JSON: `grading_order=4`, `target_reflection=1e-6`.

### DFT (probe 0 at \(x=0\))

```sh
mpirun -np 1 ./build/gnu-release-mpi/bin/opensemba_dgtd -i testData/maxwellInputs/1D_PML_L2/1D_PML_L2.json
python3 scripts/pml_dft_reflection.py Exports/single-core/1D_PML_L2 --probe 0
```

**2026-09-04:** \(R_{\mathrm{dB}}(f_{\mathrm{peak}})\approx -154~\mathrm{dB}\) → **PASS** (−40 dB gate).
See sweep notes in [`1D_PML_L3/README.md`](../1D_PML_L3/README.md).
