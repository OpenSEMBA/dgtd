# scripts/

Auxiliary post-processing tools for OpenSEMBA/dgtd (not required to build the solver).

## `pml_dft_reflection.py`

Windowed DFT reflection coefficient from PointProbe exports.

**Acceptance gate** ([docs/pml/05-verification.md](../docs/pml/05-verification.md)):

\[
20\log_{10}\frac{|E_{\mathrm{ref}}(f)|}{|E_{\mathrm{inc}}(f)|} \le -40~\mathrm{dB}
\]

at the frequency where \(|E_{\mathrm{inc}}|\) peaks in a chosen band.

### Quick start (case `1D_PML`)

Use after classical ADE-PML is implemented (CFS path is paused; see [docs/pml/27-gedney-cfs-paused.md](../docs/pml/27-gedney-cfs-paused.md)).

```sh
# from repo root, after building opensemba_dgtd
mpirun -np 1 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i testData/maxwellInputs/1D_PML/1D_PML.json

python3 scripts/pml_dft_reflection.py Exports/single-core/1D_PML \
  --probe 0 --component Ey \
  --inc-window 3.0 7.2 --ref-window 8.5 11.5 \
  --csv Exports/single-core/1D_PML/dft_P0.csv
```

Probe time in `.dat` files is SI (`t_code / c_SI`). Windows are **normalized code time** by default.

Requires: Python 3 + NumPy.
