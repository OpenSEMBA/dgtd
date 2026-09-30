# scripts/

Auxiliary post-processing tools for OpenSEMBA/dgtd (not required to build the solver).

Box PML mesh tagging and the 2D dipole absorption table: [docs/designing-a-box-pml.md](../docs/designing-a-box-pml.md).

## `pml_dft_reflection.py`

Windowed DFT reflection coefficient from PointProbe exports.

$$
20\log_{10}\frac{|E_{\mathrm{ref}}(f)|}{|E_{\mathrm{inc}}(f)|}
$$

at the frequency where $|E_{\mathrm{inc}}|$ peaks in a chosen band.

### Quick start (case `1D_PML`)

```sh
# from repo root, after building opensemba_dgtd
mpirun -np 1 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i testData/maxwellInputs/1D_PML/1D_PML.json

python3 scripts/pml_dft_reflection.py exports/SimulationData/single-core/1D_PML \
  --probe 0 --component Ey \
  --inc-window 3.0 7.2 --ref-window 8.5 11.5 \
  --csv exports/SimulationData/single-core/1D_PML/dft_P0.csv
```

Probe time in `.dat` files is SI (`t_code / c_SI`). Windows are **normalized code time** by default.

Requires: Python 3 + NumPy.

## `dipole_pml_grading_dft.py`

Same DFT, with automatic incident/reflection windows for the 2D dipole grading sweep. Replay command: [docs/designing-a-box-pml.md](../docs/designing-a-box-pml.md).

## `analyze_probe_matrix.py`

Optional check that global and `hesthaven` point-probe traces agree on paired cases. Not required to run a PML or RCS case.
