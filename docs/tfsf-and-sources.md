# Total-field / scattered-field (TFSF) and excitation sources

How `semba-dgtd` partitions the mesh into total-field (TF) and scattered-field (SF) regions, sizes and delays Gaussian pulses, and injects incident waves on `GlobalEvolution`.

---

## How TFSF splits the mesh

A Total-Field / Scattered-Field (TFSF) interface is an interior mesh boundary (a point in 1D, a closed curve in 2D, or a closed surface in 3D) whose Gmsh tag is listed under `sources[].tags` — **never** under `model.boundaries`.

At startup, `TotalFieldScatteredFieldSubMesher` (`src/components/SubMesher.cpp`) inspects the elements adjacent to the tagged TFSF faces and constructs two submeshes:

- **Total-Field (TF) side:** Contains the superposition of the incident wave and the scattered wave ($\mathbf{E}^{\mathrm{tot}} = \mathbf{E}^{\mathrm{inc}} + \mathbf{E}^{\mathrm{scat}}$). The target/scatterer sits inside the TF region.
- **Scattered-Field (SF) side:** Contains only the scattered (or reflected) wave ($\mathbf{E}^{\mathrm{scat}}$). The incident wave is canceled across the TFSF face jump so it never enters the SF region.

### 1D, 2D, and 3D TFSF topology

1. **1D TFSF (`assignIndividualTFSFAttsOnePoint1D` / `TwoPoints1D`):**
   - **Single interior point (1 point tagged):** The element to the left (`Elem1`) is tagged as **SF** and the element to the right (`Elem2`) is tagged as **TF**. A $+x$-propagating wave launches at the TFSF point, travels right through the TF region toward the target, and any reflected wave travels left across the TFSF point into the SF region (terminated by SMA or PML on the left boundary).
   - **Two points (total-field interval):** Tagging two interior points encloses a bounded TF interval between them, with SF outside on both ends.
2. **2D and 3D closed TFSF surface:**
   - The tagged curves (2D) or surface triangles/quads (3D) form a closed Huygens box around the scatterer.
   - `SubMesher` oriented-normal checks tag the interior elements adjacent to the box as **TF** and the exterior elements as **SF**.
   - In MPI runs, the partitioner (`applyMetisPartitioningWithTFSFComponents` in `src/driver/driver.cpp`) locks each pair of elements sharing a TFSF face onto the same MPI rank so the TFSF jump operator is strictly rank-local.

### Where to place probes, targets, and absorbers

```text
+-------------------------------------------------------------------+
| Outer boundary: SMA or volumetric Box PML                         |
|                                                                   |
|   Scattered-Field (SF) region                                     |
|   - Place reflected-only PointProbes here (1D)                    |
|   - Place closed rcssurface / farfield Huygens surface here       |
|                                                                   |
|   +-----------------------------------------------------------+   |
|   | TFSF interface (tags listed in sources[].tags)            |   |
|   |                                                           |   |
|   |   Total-Field (TF) region                                 |   |
|   |   - Place scatterer (PEC, SGBC, dielectric, plasma) here  |   |
|   |   - Place incoming/total-field PointProbes here           |   |
|   +-----------------------------------------------------------+   |
+-------------------------------------------------------------------+
```

---

## Gaussian pulse width (`f_1e` vs. `spread`) and auto-delay

Time-domain Gaussian envelopes have the form:

$$
g(t') = \exp\!\left(-\frac{(t' - t_0')^2}{2\sigma^2}\right),
$$

where $t' = c_0 t$ and $\sigma$ (`spread`) are in light-metres.

### Setting bandwidth (`f_1e` or `spread`)

You can specify the Gaussian width using **either**:

- **`f_1e` (Hz):** The frequency where the incident **power** spectrum $|G(\omega)|^2 / |G(0)|^2 = \exp(-(\omega\sigma)^2)$ drops to $1/e$ ($\approx -4.34\,\mathrm{dB}$) of its DC value:

$$
\sigma = \frac{c_0}{2\pi f_{1e}}.
$$

- **`spread` (light-metres):** Sets $\sigma$ directly in metres.

At $f = 2 f_{1e}$, the spectral amplitude is $e^{-2} \approx 0.135$ ($-17.4\,\mathrm{dB}$); at $f = 3 f_{1e}$, it is $e^{-4.5} \approx 0.011$ ($-39.1\,\mathrm{dB}$). If both `f_1e` and `spread` are present in the same JSON object, `f_1e` takes precedence and a warning is printed.

### Automatic arrival delay (`AUTO_DELAY_N_SIGMA = 5.0`)

When `mean` is omitted on a `planewave` or `dipole` source, `driver.cpp` scans all boundary elements on the TFSF surface to find the earliest arrival phase along the propagation direction $\hat{\mathbf{d}}$:

$$
\phi_{\min} = \min_{\mathbf{x} \in \Gamma_{\mathrm{TFSF}}} (\mathbf{x}\cdot\hat{\mathbf{d}}), \qquad \mu_{\mathrm{1D}} = \phi_{\min} - 5\sqrt{2}\,\sigma.
$$

This places the pulse peak $5\sqrt{2}\,\sigma$ upstream of the earliest TFSF point at $t = 0$, guaranteeing that the initial field on the TFSF surface is $\exp(-25) \approx 1.4 \times 10^{-11}$ (machine-zero turn-on transient).

### Startup mesh resolution check

At solver initialization, `warnGaussianPulseVsMesh` (`src/solver/Solver.cpp`) compares the source's effective $1/e$ frequency $f_{1e} + f_{\mathrm{carrier}}$ against the mean element size $h_{\mathrm{avg}}$ and polynomial order $p$ using the rule of thumb:

$$
N_{\mathrm{ppw}} = \max\!\left(3,\, \frac{15}{p + 1}\right)\text{ elements per wavelength}, \qquad f_{\mathrm{mesh}} = \frac{c_0}{N_{\mathrm{ppw}}\,h_{\mathrm{avg}}}.
$$

If $f_{1e} + f_{\mathrm{carrier}} > f_{\mathrm{mesh}}$, the solver prints a non-fatal warning suggesting a lower `f_1e`, a larger `spread`, or mesh refinement.

---

## Supported source types

All sources are listed in the top-level `sources` array; multiple sources superimpose linearly.

### 1. Plane wave (`"type": "planewave"`)

Injects a uniform plane wave across the TFSF surface `tags` with electric polarization $\mathbf{p}$ (`polarization`), propagation direction $\hat{\mathbf{d}}$ (`propagation`), and magnetic field $\mathbf{H}^{\mathrm{inc}} = \hat{\mathbf{d}} \times \mathbf{E}^{\mathrm{inc}}$ (normalized $Z_0 = 1$).

Three excitation waveforms are supported under `magnitude`:

1. **Unmodulated Gaussian (`f_1e` or `spread`, no `frequency`):** Broadband baseband pulse from DC up to $\sim 2.5 f_{1e}$.
2. **Modulated Gaussian (`f_1e` or `spread`, plus `frequency` in Hz):** Multiplies the Gaussian envelope by $\cos\big(2\pi (f_0 / c_0)(\mathbf{x}\cdot\hat{\mathbf{d}} - t' - \mu_{\mathrm{1D}})\big)$, centering the spectral content around carrier $f_0$.
3. **Sampled excitation file (`file`):** Reads a two-column ASCII file (`time_seconds amplitude`) and interpolates linearly. Time zero corresponds to the plane $\mathbf{x}\cdot\hat{\mathbf{d}} = 0$. Do not combine `file` with `spread`, `f_1e`, `frequency`, or `mean`.

```json
{
  "type": "planewave",
  "tags": [3],
  "polarization": [0.0, 1.0, 0.0],
  "propagation": [1.0, 0.0, 0.0],
  "magnitude": {
    "type": "gaussian",
    "f_1e": 2.0e9,
    "frequency": 6.0e9
  }
}
```

### 2. Hertzian dipole (`"type": "dipole"`)

Injects the exact near-, intermediate-, and far-field radiation of an infinitesimal electric dipole of geometric scale `length` driven by a derivative-Gaussian current pulse across the closed TFSF surface `tags`.

| Field | Default | Description |
|-------|---------|-------------|
| `tags` | — | TFSF interface tags. |
| `magnitude.length` | — | Dipole length $dl$. |
| `magnitude.spread` / `f_1e` | — | Gaussian width $\sigma$ (light-metres) or $1/e$ power frequency (Hz). |
| `magnitude.amplitude_peak` | `1.0` | Target peak equatorial $|E_\theta|$ at distance `peak_radius`. |
| `magnitude.peak_radius` | auto | Reference radius for `amplitude_peak` (defaults to minimum radius on the TFSF tags). |
| `magnitude.mean` | auto | Optional pulse delay override. |

### 3. Delta-gap slot source (`"type": "delta_gap"`)

Applies an impressed tangential electric field across an interior face (`tags`: a curve in 2D or a rectangular surface in 3D) with vacuum on both sides. Unlike `planewave` or `dipole`, `delta_gap` uses its own decoupled same-side Zero/Two-normal face matrix rather than a TFSF jump.

- Both sides of the gap face receive $+\tfrac{1}{2}\mathbf{E}_{\mathrm{pol}}\,w(t')$ and $\mathbf{H}^{\mathrm{inc}} = \mathbf{0}$.
- `signal` selects `"gaussian"` (default) or `"gaussian_derivative"`, centered at $t_0' = 5\sigma$.
- If both `spread` and `f_1e` are omitted, $\sigma$ is chosen automatically from the gap length $L$ so the spectrum is `db_cut` (default `-20` dB) down at $\mathrm{magnitude} \times L / 10$.
- In 3D, the solver measures the slot separation $h$ (along `polarization`) and width $w$ (perpendicular in the face) and prints the wide-plate line impedance $Z = h/w$ and capacitance per unit length $C' = w/h$.

### 4. Coaxial TEM port (`"type": "coaxial_port"`)

Launches a circular coaxial TEM mode across an annular port face (`tags.load`) bounded by inner and outer PEC conductors (`tags.live` and `tags.outer`) in 3D:

$$
E_\rho(\rho, t') = \frac{V_0\,w(t')}{\rho\ln(b/a)}, \qquad \mathbf{H} = \hat{\mathbf{s}} \times \mathbf{E},
$$

where $a$ and $b$ are the inner and outer conductor radii measured automatically from the mesh vertices on `tags.live` and `tags.outer`, $V_0$ is `magnitude` (default `1.0`), and $\hat{\mathbf{s}}$ points from the SF volume into the TF volume.

- **Interior annular load face:** Acts as a TFSF interface. The solver automatically identifies the SF neighbor (the side touching an SMA boundary or closer to a PML volume) and the TF neighbor (the coaxial transmission line feeding the device).
- **Boundary annular load face:** The `tags.load` tag must also be listed as `"SMA"` in `model.boundaries`; the face acts as a matched TEM port launching into the single coaxial volume and absorbing outgoing waves.

### 5. Volumetric initial condition (`"type": "initial"`)

Projects an analytical spatial distribution onto $\mathbf{E}$ or $\mathbf{H}$ at $t = 0$ (`SourcesManager::setInitialFields`), used for cavity resonance and convergence benchmarks:

- `"gaussian"`: spatial Gaussian $\exp(-|\mathbf{x} - \mathbf{x}_0|^2 / (2\sigma^2))$ over `dimension` ($1, 2,$ or $3$) active coordinates.
- `"resonant"`: standing-wave sinusoidal cavity mode specified by `magnitude.modes`.
- `"besselj6_2D"` / `"besselj6_3D"`: cylindrical and spherical $J_6$ Bessel cavity eigenmodes.
