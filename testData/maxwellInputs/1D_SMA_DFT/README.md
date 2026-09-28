# 1D_SMA_DFT

SMA-only baseline comparable to [`1D_PML`](../1D_PML/): materials are **all vacuum**
(no PML). Outer BCs remain SMA.

Use with [`scripts/pml_dft_reflection.py`](../../../scripts/pml_dft_reflection.py) to
compare reflection against the volumetric PML case. Box PML tagging and the 2D
absorption table: [`docs/designing-a-box-pml.md`](../../../docs/designing-a-box-pml.md).
