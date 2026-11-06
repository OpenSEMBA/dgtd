# RCS and near-to-far-field (NTFF) post-processing

How to record scattered surface fields during a 2D or 3D `opensemba_dgtd` simulation and compute bistatic or monostatic Radar Cross Section (RCS) and far-field radiation patterns with `opensemba_rcs`.

---

## Two-stage RCS workflow

`semba-dgtd` decouples time-domain surface recording from frequency- and angle-domain far-field integration:

1. **Stage 1 — Time-domain surface dump (`opensemba_dgtd`):**
   Add a `probes.rcssurface` entry to the Maxwell case JSON. During the run, `RCSSurfaceExporter` extracts the tangential $(\mathbf{E}, \mathbf{H})$ fields on the tagged Huygens surface and writes a compact binary dump (`surface_data.bin`) per MPI rank:

```text
exports/SimulationData/<runmode>/<casename>/RCSSurface/<probe_name>/
└── rank0/
    ├── mesh
    └── surface_data.bin
```

2. **Stage 2 — Offline DFT and surface integration (`opensemba_rcs`):**
   Point `opensemba_rcs` at a `*.rcs.json` sweep file. You can re-run `opensemba_rcs` as many times as needed with different frequency bands, angular resolutions, time windows (`max_time`), or subsampling strides (`every_n_steps`) without re-running the time-domain solver:

```sh
mpiexec -n 1 ./build/gnu-release-mpi/bin/opensemba_rcs \
  -i testData/rcsInputs/3D_RCS_SGBC_CoatedSphere_O2.rcs.json
```

### `probes.rcssurface` vs. `probes.farfield`

- **`probes.rcssurface` + `opensemba_rcs` (Recommended):** Stores binary surface snapshots once (`surface_data.bin`), supports MPI-distributed surfaces (`rank0/`, `rank1/`, …), automatically switches to RAM-gated streaming DFT or CUDA GPU streaming DFT for large 3D runs, and normalizes RCS by the analytical incident plane-wave spectrum.
- **`probes.farfield`:** Legacy online NTFF probe (`NearToFarFieldProbes/<name>/`) that assembles `ParLinearForm` face integrals (`NTFFBdrFaceIntegrator`) directly from the case JSON. New 2D/3D scattering benchmarks should use `probes.rcssurface`.

---

## Where to place the Huygens surface in the mesh

The boundary tags listed in `probes.rcssurface[].tags` must form a **closed interior surface** (a closed loop of curves in 2D, or a closed shell of triangles/quads in 3D) enclosing the target:

```text
Outer Absorber (SMA or Box PML)
  └── Scattered-Field (SF) region
        └── Closed rcssurface Huygens shell (probes.rcssurface[].tags)
              └── Closed TFSF excitation surface (sources[].tags)
                    └── Scatterer (PEC, SGBC, dielectric, Lorentz/Drude)
```

- **Standard placement (inside SF):** Place the `rcssurface` shell in the Scattered-Field region, outside the TFSF surface and inside the PML/SMA absorber. Outward normals $\hat{\mathbf{n}}$ are oriented automatically by `NearToFarFieldSubMesher` (`src/components/SubMesher.cpp`).
- **Temporal sampling (`steps` or `saves`):** Set `steps` (e.g. `"steps": 1` or `"steps": 5`) so that the sampling interval $\Delta t_{\mathrm{export}} = \mathrm{steps} \times \Delta t / c_0$ provides at least $10\text{--}20$ samples per period at the highest frequency in your RCS sweep.

Case JSON example (`testData/maxwellInputs/3D_RCS_SGBC_CoatedSphere_O2/3D_RCS_SGBC_CoatedSphere_O2.json`):

```json
"probes": {
  "rcssurface": [
    {
      "name": "RCSSurfaceProbe",
      "tags": [3],
      "steps": 1
    }
  ]
}
```

---

## How `RCSSurfacePostProcessor` computes far-field and RCS

For each requested frequency $f$ ($k = 2\pi f / c_0$) and observation direction $\hat{\mathbf{r}}(\theta, \phi)$:

1. **Discrete Fourier Transform (DFT):**
   Transforms the time-domain surface fields $(\mathbf{E}(\mathbf{r}', t), \mathbf{H}(\mathbf{r}', t))$ at every surface quadrature point $\mathbf{r}'$ into phasors $(\tilde{\mathbf{E}}(\mathbf{r}', f), \tilde{\mathbf{H}}(\mathbf{r}', f))$ using trapezoidal time weights.
2. **Equivalent surface currents:**
   Computes the electric and magnetic Love equivalence currents on the Huygens surface:

$$
\mathbf{J}_s(\mathbf{r}', f) = \hat{\mathbf{n}}(\mathbf{r}') \times \tilde{\mathbf{H}}(\mathbf{r}', f), \qquad \mathbf{M}_s(\mathbf{r}', f) = -\hat{\mathbf{n}}(\mathbf{r}') \times \tilde{\mathbf{E}}(\mathbf{r}', f).
$$

3. **Far-field radiation potentials:**
   Integrates the surface currents over all quadrature points across all MPI rank dumps with the far-field phase factor $\exp(+jk\,\mathbf{r}'\cdot\hat{\mathbf{r}})$ to obtain the transverse far-field components $(E_\theta^{\mathrm{ff}}, E_\phi^{\mathrm{ff}})$.
4. **Incident spectrum normalization:**
   Reads the `planewave` source from `testData/maxwellInputs/<casename>/<casename>.json` (`Gaussian` or `ModulatedGaussian`), evaluates its analytical Fourier transform $|\tilde{E}_{\mathrm{inc}}(f)|^2$, and normalizes the radar cross section:
   - **3D RCS ($\mathrm{m}^2$):**

$$
\sigma_{\mathrm{3D}}(\theta, \phi, f) = 4\pi\,\frac{|E_\theta^{\mathrm{ff}}|^2 + |E_\phi^{\mathrm{ff}}|^2}{|\tilde{E}_{\mathrm{inc}}(f)|^2}.
$$

   - **2D scattering width ($\mathrm{m}$):**

$$
\sigma_{\mathrm{2D}}(\phi, f) = 2\pi\,\frac{|E_\theta^{\mathrm{ff}}|^2 + |E_\phi^{\mathrm{ff}}|^2}{k\,|\tilde{E}_{\mathrm{inc}}(f)|^2}.
$$

### Spherical angle conventions

Angles in `*.rcs.json` are specified in **radians** (`{ "start": ..., "end": ..., "steps": ... }`):

- **$\theta \in [0, \pi]$:** Polar angle measured from the $+Z$ axis ($\theta = 0$ is $+Z$, $\theta = \pi/2$ is the $XY$ equatorial plane, $\theta = \pi$ is $-Z$).
- **$\phi \in [0, 2\pi)$:** Azimuthal angle in the $XY$ plane measured from $+X$ toward $+Y$:

$$
\hat{\mathbf{r}}(\theta, \phi) = (\sin\theta\cos\phi,\; \sin\theta\sin\phi,\; \cos\theta).
$$

- **2D restriction:** In 2D simulations, scattering occurs in the $XY$ plane, so `angles.theta` **must** be fixed at $\theta = \pi/2$ (`1.5707963267948966` with `"steps": 1`).

---

## RAM gating and CUDA streaming DFT

Large 3D scattering runs (such as the NASA Almond at $5\,\mathrm{GHz}$) can produce tens of gigabytes of surface snapshots across multiple MPI ranks. `RCSSurfacePostProcessor` automatically manages memory per rank directory:

1. **RAM budget estimation:**
   Before loading `surface_data.bin`, `peekSurfaceFile` reads the geometry header and counts the kept snapshots $N_{\mathrm{kept}}$. It estimates the peak memory required to hold all snapshots in RAM (`estimateLoadAllBytes`) and compares it against the budget:
   - If `"ram_gate"` (in GiB) is set in `*.rcs.json`, that value is used as the hard RAM budget.
   - If `"ram_gate"` is omitted, the budget defaults to **$50\%$ of Linux `MemAvailable`** (`/proc/meminfo`).
2. **Load-all vs. streaming DFT:**
   - **Under budget (load-all):** Reads all kept snapshots into memory once and computes the DFT.
   - **Over budget (streaming):** Streams `surface_data.bin` snapshot-by-snapshot into the frequency-domain accumulator (`dftStreamFromFile`), keeping peak memory proportional to $N_{\mathrm{freq}} \times N_{\mathrm{dofs}}$ regardless of the number of time steps.
3. **CUDA GPU streaming DFT (`tryCudaStreamDft`):**
   When `opensemba_rcs` is built with CUDA (`gnu-release-cuda`), streaming DFT automatically runs on the GPU if the complex frequency accumulator fits in free VRAM, falling back to the CPU streaming path otherwise.

---

## Offline `opensemba_rcs` JSON example and output files

Example (`testData/rcsInputs/3D_RCS_SGBC_CoatedSphere_O2.rcs.json`):

```json
{
  "runmode": "cuda-1",
  "casename": "3D_RCS_SGBC_CoatedSphere_O2",
  "frequencies": {
    "start": 100000000.0,
    "end": 1500000000.0,
    "steps": 141
  },
  "angles": {
    "theta": {
      "start": 3.141592653589793,
      "end": 3.141592653589793,
      "steps": 1
    },
    "phi": {
      "start": 0.0,
      "end": 0.0,
      "steps": 1
    }
  }
}
```

Outputs are written under `exports/SimulationData/<runmode>/<casename>/RCSSurface/<probe_name>/`:

- `farfield/farfieldData_Th_<theta>_Phi_<phi>_dgtd.dat` — unnormalized far-field radiation potential vs. frequency (Hz).
- `rcs/rcsData_Th_<theta>_Phi_<phi>_dgtd.dat` — normalized RCS ($\mathrm{m}^2$ in 3D, $\mathrm{m}$ in 2D) vs. frequency (Hz).

For MOR-based RCS post-processing from reduced-order bases (`Ur.bin` and `xr/`), see [`docs/mor2rcs.md`](./mor2rcs.md).
