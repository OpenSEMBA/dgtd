# 1D_SMA_DFT

SMA-only baseline comparable to [`1D_PML`](../1D_PML/): materials are **all vacuum**
(no PML). Outer BCs remain SMA.

Use with [`scripts/pml_dft_reflection.py`](../../../scripts/pml_dft_reflection.py) to
compare reflection against the volumetric PML case once classical ADE-PML is active
(see [`docs/pml/27-gedney-cfs-paused.md`](../../../docs/pml/27-gedney-cfs-paused.md)).

### Historical CFS-era note

During the paused Gedney CFS attempt, SMA-alone on this mesh was ≪ −100 dB while the
CFS PML case still showed interface reflection above −40 dB on some windows. Re-baseline
after classical ADE lands.
