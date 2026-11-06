# Documentation index

| Document | Description |
|----------|-------------|
| [json-input-format.md](./json-input-format.md) | Full JSON schema reference: Maxwell case files (materials, Debye/Lorentz/Drude, spatial profiles, boundaries, probes, sources) and offline `opensemba_rcs` inputs |
| [designing-a-box-pml.md](./designing-a-box-pml.md) | How to cut, tag, and grade a Cartesian box SC-PML, with the 2D absorption reference |
| [dispersive-materials.md](./dispersive-materials.md) | Single-pole Debye, Lorentz, and cold-plasma Drude formulations, spatial rate profiles (`gaussian`, `positive_gaussian_derivative`, `ramc_61km_*`), and 1D automatic region stretching |
| [sgbc.md](./sgbc.md) | Surface general boundary condition (SGBC): multi-layer thin sheets and RAM coatings, automatic layer meshing, IMEX sub-stepping, and ParaView export |
| [tfsf-and-sources.md](./tfsf-and-sources.md) | Total-field / scattered-field (TFSF) mesh partitioning, Gaussian `f_1e` / `spread` and auto-delay, and `planewave`, `dipole`, `delta_gap`, `coaxial_port`, and `initial` excitations |
| [rcs-and-farfield.md](./rcs-and-farfield.md) | Two-stage RCS and near-to-far-field workflow (`probes.rcssurface` + `opensemba_rcs`), Huygens surface placement, angle conventions, and RAM-gated/CUDA streaming DFT |
| [mor2paraview.md](./mor2paraview.md) | MOR state $\to$ ParaView replay tool (`opensemba_mor2paraview`) |
| [mor2rcs.md](./mor2rcs.md) | MOR reconstruct $\to$ RCS post-processing (`opensemba_mor2rcs`) |
| [mfem-usage.md](./mfem-usage.md) | How the solver uses MFEM: stock library objects, and the objects this solver specializes (TFSF submeshes, integrators, fork hooks) |
