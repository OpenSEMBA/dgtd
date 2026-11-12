# Dispersive materials and spatial profiles

Volumetric single-pole electric Debye, Lorentz, and cold-plasma Drude materials on `GlobalEvolution` and `HesthavenEvolution`.

Internal solver units are normalized ($c_0 = \varepsilon_0 = \mu_0 = 1$, time in light-metres $t' = c_0 t$). JSON rates and relaxation times are given in SI units (seconds, Hz, or rad/s) and converted at parse time. Use `"evolution_operator": "global"` (or omit it) or `"hesthaven"`. `"ode_type": 0` is RK4. Implicit ODE types include these poles. `"hesthaven"` is 2D/3D and applies the same ADE as element-local updates.

---

## Material formulations

### 1. Single-pole electric Debye (`debye`)

Models polar dielectrics with frequency-dependent relative permittivity:

$$
\varepsilon_r(\omega) = \varepsilon_\infty + \frac{\varepsilon_s - \varepsilon_\infty}{1 + j\omega\tau}.
$$

In the time domain, `GlobalEvolution` appends $3N$ polarization degrees of freedom $(P_x, P_y, P_z)$ and integrates the auxiliary differential equation (ADE) alongside Maxwell's curl equations:

$$
\tau'\frac{\partial \mathbf{P}}{\partial t'} + \mathbf{P} = (\varepsilon_s - \varepsilon_\infty)\,\mathbf{E}, \qquad \varepsilon_\infty\frac{\partial \mathbf{E}}{\partial t'} = \nabla \times \mathbf{H} - \sigma_{\mathrm{bulk}}\mathbf{E} - \frac{\partial \mathbf{P}}{\partial t'},
$$

where $\tau' = c_0 \tau_{\mathrm{SI}}$ is the relaxation time in light-metres.

| Field | Type | Requirement | Description |
|-------|------|-------------|-------------|
| `eps_inf` | double | $\ge 1.0$ | Instantaneous (high-frequency) relative permittivity $\varepsilon_\infty$. |
| `eps_s` | double | $> \text{eps\_inf}$ | Static (DC) relative permittivity $\varepsilon_s$. |
| `tau` | double | $> 0$ | Relaxation time $\tau_{\mathrm{SI}}$ in seconds. |

Example (`testData/maxwellInputs/1D_Debye/1D_Debye.json`):

```json
{
  "tags": [3],
  "debye": {
    "eps_inf": 2.0,
    "eps_s": 5.0,
    "tau": 5.308837458876128e-10
  }
}
```

### 2. Single-pole electric Lorentz (`lorentz`)

Models resonant dispersive media and magnetoplasmas/cold plasmas with relative permittivity:

$$
\varepsilon_r(\omega) = \varepsilon_\infty + \frac{\omega_p^2}{\omega_1^2 - \omega^2 + j\gamma\omega}.
$$

In the time domain, `GlobalEvolution` appends $6N$ auxiliary degrees of freedom $(P_x, P_y, P_z, J_x, J_y, J_z)$ and integrates the second-order pole as a first-order system:

$$
\frac{\partial \mathbf{P}}{\partial t'} = \mathbf{J}, \qquad \frac{\partial \mathbf{J}}{\partial t'} = \tilde{\omega}_p^2(\mathbf{x})\,\mathbf{E} - \tilde{\omega}_1^2(\mathbf{x})\,\mathbf{P} - \tilde{\gamma}(\mathbf{x})\,\mathbf{J}, \qquad \varepsilon_\infty\frac{\partial \mathbf{E}}{\partial t'} = \nabla \times \mathbf{H} - \sigma_{\mathrm{bulk}}\mathbf{E} - \mathbf{J},
$$

where $\tilde{\omega}_p = \omega_p / c_0$, $\tilde{\omega}_1 = \omega_1 / c_0$, and $\tilde{\gamma} = \gamma / c_0$ have units of $\mathrm{m}^{-1}$.

Rates can be specified in **either** ordinary frequencies in Hz (`f_p`, `f_1`, `gamma` / `nu_c`) **or** angular frequencies in rad/s (`omega_p`, `omega_1`, `gamma` / `nu_c`), never mixed in the same block:

- **Hz style:** `f_p`, `f_1`, and `gamma` (or `nu_c`) are all multiplied by $2\pi$ before dividing by $c_0$.
- **rad/s style:** `omega_p`, `omega_1`, and `gamma` (or `nu_c`) are divided directly by $c_0$.

| Field | Type | Requirement | Description |
|-------|------|-------------|-------------|
| `eps_inf` | double | $\ge 1.0$ | Instantaneous relative permittivity $\varepsilon_\infty$. |
| `f_p` / `omega_p` | double or object | $> 0$ | Plasma frequency in Hz (`f_p`) or rad/s (`omega_p`), or spatial profile object. |
| `f_1` / `omega_1` | double or object | $\ge 0$ | Resonance frequency in Hz (`f_1`) or rad/s (`omega_1`), or spatial profile object ($\omega_1 = 0$ is a cold plasma). |
| `gamma` / `nu_c` | double or object | $\ge 0$ | Collision/damping rate in Hz (with `f_*`) or rad/s / $\mathrm{s}^{-1}$ (with `omega_*`), or spatial profile object. |

### 3. Cold-plasma Drude (`drude`)

`drude` is a dedicated shorthand for a cold-plasma Lorentz pole with $\omega_1 \equiv 0$:

$$
\varepsilon_r(\omega) = \varepsilon_\infty - \frac{\omega_p^2}{\omega^2 - j\nu_c\omega}.
$$

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `eps_inf` | double | `1.0` | Instantaneous relative permittivity ($\ge 1.0$). |
| `omega_p` / `f_p` | double or object | — | Plasma frequency in rad/s (`omega_p`) or Hz (`f_p`), or spatial profile object. |
| `nu_c` / `gamma` | double or object | — | Collision frequency in $\mathrm{s}^{-1}$ (with `omega_p`) or Hz (with `f_p`), or spatial profile object. |

### Material block rules

- Do **not** set `relative_permittivity` on a `debye`, `lorentz`, or `drude` material; the electric mass matrix $M_E$ automatically uses `eps_inf`.
- `relative_permeability` (default `1.0`) and `bulk_conductivity` (default `0.0` S/m) remain available alongside the dispersive pole.
- A mesh tag cannot combine `debye`, `lorentz`, `drude`, or `type: "PML"`, and SGBC layers cannot carry dispersive poles.

---

## Spatial rate profiles and DG mass assembly

### Weak-form quadrature integration

When all rates on an element $K$ are constant scalars, the element mass matrix cancels against $M_K^{-1}$ and the auxiliary coupling block reduces to a point-diagonal scaling on each degree of freedom.

When any rate ($\omega_p$, $\omega_1$, or $\nu_c$) is given as a spatial profile object, the operator wraps the profile in an `mfem::FunctionCoefficient`, integrates it at the element's Gauss quadrature points via `mfem::MassIntegrator`, and pre-multiplies by the element inverse mass matrix $(M_K^{(1)})^{-1}$:

$$
\frac{\mathrm{d}\mathbf{J}_K}{\mathrm{d}t'} = (M_K^{(1)})^{-1} M_K^{(\tilde{\omega}_p^2)}\,\mathbf{E}_K - (M_K^{(1)})^{-1} M_K^{(\tilde{\omega}_1^2)}\,\mathbf{P}_K - (M_K^{(1)})^{-1} M_K^{(\tilde{\gamma})}\,\mathbf{J}_K,
$$

where

$$
\big(M_K^{(f)}\big)_{ij} = \int_K f(\mathbf{x})\,\phi_i(\mathbf{x})\,\phi_j(\mathbf{x})\,\mathrm{d}\mathbf{x}.
$$

Evaluating $\tilde{\omega}_p^2(\mathbf{x}) \ge 0$ and $\tilde{\gamma}(\mathbf{x}) \ge 0$ directly at quadrature points avoids polynomial interpolation undershoot (Runge/Gibbs oscillations) across steep shock fronts and preserves element-block-diagonal sparsity with zero MPI halo overhead.

### 1D automatic region stretching

In 1D simulations, spatial profiles are not locked to hardcoded mesh coordinates. At initialization, the solver scans the mesh for all elements assigned to the material's `tags` and records the physical bounding interval $[x_{\min}^{\mathrm{mesh}},\, x_{\max}^{\mathrm{mesh}}]$. Every quadrature point $x \in [x_{\min}^{\mathrm{mesh}},\, x_{\max}^{\mathrm{mesh}}]$ is mapped to the normalized unit coordinate:

$$
\xi(x) = \frac{x - x_{\min}^{\mathrm{mesh}}}{x_{\max}^{\mathrm{mesh}} - x_{\min}^{\mathrm{mesh}}} \in [0, 1].
$$

Each profile evaluates its canonical shape on $\xi \in [0, 1]$, so the complete profile automatically stretches or compresses to cover the tagged mesh region whether the region is $10\,\mathrm{mm}$, $12\,\mathrm{mm}$, $24\,\mathrm{mm}$, or $2\,\mathrm{m}$ wide.

### Built-in profile types

| `type` | Required keys | Optional keys (defaults) | Shape on $\xi \in [0, 1]$ |
|--------|---------------|--------------------------|---------------------------|
| `"constant"` | `value` | — | Uniform value across the region. |
| `"gaussian"` | `peak` (or `amplitude`) | `cutoff_factor` (`1e-2`), `offset` (`0.0`) | Symmetric bell curve from $\epsilon\cdot\mathrm{peak}$ at $\xi=0$ to $\mathrm{peak}$ at $\xi=0.5$ and back to $\epsilon\cdot\mathrm{peak}$ at $\xi=1$. |
| `"positive_gaussian_derivative"` | `peak` (or `amplitude`) | `cutoff_factor` (`1e-2`), `offset` (`0.0`) | Two-lobe rectified Gaussian derivative with unit peaks at $z = \pm 1$, central null at $\xi=0.5$, and outer tails $\epsilon\cdot\mathrm{peak}$ at $\xi \in \{0, 1\}$. |
| `"ramc_61km_omega_p"` | — | `direction` (`1.0`) | RAM-C II 61 km stagnation-line $\omega_p$ $[\mathrm{rad/s}]$ stretched from $[-0.012, 0.0]\,\mathrm{m}$ onto $\xi \in [0, 1]$. |
| `"ramc_61km_nu_c_gupta"` | — | `direction` (`1.0`) | RAM-C II 61 km stagnation-line $\nu_c$ $[\mathrm{s}^{-1}]$ (Gupta closure) stretched onto $\xi \in [0, 1]$. |
| `"ramc_61km_nu_c_liu"` | — | `direction` (`1.0`) | RAM-C II 61 km stagnation-line $\nu_c$ $[\mathrm{s}^{-1}]$ (Liu closure) stretched onto $\xi \in [0, 1]$. |
| `"ramc_50km_omega_p"` | — | `direction` (`1.0`) | RAM-C II 50 km stagnation-line $\omega_p$ $[\mathrm{rad/s}]$ stretched from $[-0.024, 0.0]\,\mathrm{m}$ onto $\xi \in [0, 1]$. |
| `"ramc_50km_nu_c_gupta"` | — | `direction` (`1.0`) | RAM-C II 50 km stagnation-line $\nu_c$ $[\mathrm{s}^{-1}]$ (Gupta closure) stretched onto $\xi \in [0, 1]$. |
| `"ramc_50km_nu_c_liu"` | — | `direction` (`1.0`) | RAM-C II 50 km stagnation-line $\nu_c$ $[\mathrm{s}^{-1}]$ (Liu closure) stretched onto $\xi \in [0, 1]$. |

#### Mathematical definition of `"gaussian"` and `"positive_gaussian_derivative"`

Given `cutoff_factor` $\epsilon \in (0, 1)$ (default `1e-2`), both profiles determine the canonical cutoff coordinate $z_{\mathrm{cut}}$ where the tail decays to $\epsilon$ of the peak, and map $[-z_{\mathrm{cut}}, +z_{\mathrm{cut}}]$ onto $\xi \in [0, 1]$ via $z(\xi) = (2\xi - 1)\,z_{\mathrm{cut}}$:

1. **`"gaussian"`:**

$$
z_{\mathrm{cut}} = \sqrt{-2\ln(\epsilon)} \quad (\approx 3.035\text{ for }\epsilon = 10^{-2}),
$$

$$
f(\xi) = \mathrm{offset} + \mathrm{peak}\cdot\exp\!\left(-\frac{z(\xi)^2}{2}\right).
$$

2. **`"positive_gaussian_derivative"`:**

Uses the non-negative normalized derivative envelope $D(z) = |z|\,\exp\!\big(\tfrac{1}{2}(1 - z^2)\big)$, which has two symmetric lobes peaking at $D(\pm 1) = 1$ and a central trough $D(0) = 0$. The outer cutoff $z_{\mathrm{cut}} > 1$ is solved via Newton–Raphson from:

$$
z_{\mathrm{cut}}\,\exp\!\left(\frac{1 - z_{\mathrm{cut}}^2}{2}\right) = \epsilon \quad (\approx 3.648\text{ for }\epsilon = 10^{-2}),
$$

$$
f(\xi) = \mathrm{offset} + \mathrm{peak}\cdot |z(\xi)|\,\exp\!\left(\frac{1 - z(\xi)^2}{2}\right).
$$

Across $\xi \in [0, 1]$:
- $\xi = 0$ and $\xi = 1$: outer tail $\mathrm{offset} + \epsilon\cdot\mathrm{peak}$
- $\xi = \tfrac{1}{2}(1 - 1/z_{\mathrm{cut}}) \approx 0.363$: first lobe peak $\mathrm{offset} + \mathrm{peak}$
- $\xi = 0.5$: central trough $\mathrm{offset}$
- $\xi = \tfrac{1}{2}(1 + 1/z_{\mathrm{cut}}) \approx 0.637$: second lobe peak $\mathrm{offset} + \mathrm{peak}$

#### RAM-C II 61 km stagnation-line profiles

The presets `"ramc_61km_omega_p"`, `"ramc_61km_nu_c_gupta"`, and `"ramc_61km_nu_c_liu"` evaluate the analytical curve fits for the RAM-C II flight experiment at 61 km altitude along the $12\,\mathrm{mm}$ nose stagnation line ($s \in [-0.012, 0.0]\,\mathrm{m}$, with the bow shock near $s \approx -7.31\,\mathrm{mm}$, peak plasma frequency $f_{p,\max} \approx 108.1\,\mathrm{GHz}$, and vehicle nose at $s = 0$):

$$
f(s) = f_\infty + \frac{1}{2}\left(1 + \tanh\!\left(\frac{s - x_{\mathrm{sh}}}{d_{\mathrm{sh}}}\right)\right)\max\!\left(0,\, a_0 + a_1\tilde{s} + a_2\tilde{s}^2 + a_3\tilde{s}^3 \pm A_{\mathrm{bl}}\exp\!\left(\frac{s}{d_{\mathrm{bl}}}\right)\right),
$$

with $\tilde{s} = s + 0.0066\,\mathrm{m}$ and $s(\xi) = -0.012\,(1 - \xi)$ for `direction: 1.0` (or $s(\xi) = -0.012\,\xi$ for `direction: -1.0` if the nose sits at the left end of the region).

The 50 km presets `"ramc_50km_omega_p"`, `"ramc_50km_nu_c_gupta"`, and `"ramc_50km_nu_c_liu"` use the same shock-plus-polynomial shape on the $24\,\mathrm{mm}$ stagnation line $s \in [-0.024, 0.0]\,\mathrm{m}$ ($s(\xi) = -0.024\,(1 - \xi)$ when `direction` is $1$). Peak $f_p$ on that fit is about $0.96\,\mathrm{GHz}$ at the nose. The 1D check is [`testData/testInputs/1D_Drude_RAMC_50km/1D_Drude_RAMC_50km.json`](../testData/testInputs/1D_Drude_RAMC_50km/1D_Drude_RAMC_50km.json) with [`scripts/ramc_50km_s11_skrf.py`](../scripts/ramc_50km_s11_skrf.py), from $1\,\mathrm{MHz}$ to $5\,\mathrm{GHz}$.

---

## Time-step and mesh resolution guidelines

1. **Plasma oscillation time-step constraint:**
   In explicit RK4 integration, the plasma oscillation frequency $\tilde{\omega}_{p,\max} = \omega_{p,\max} / c_0$ must satisfy the RK4 imaginary-axis stability limit $\Delta t' \cdot \tilde{\omega}_{p,\max} \le 2\sqrt{2} \approx 2.83$, in addition to the standard DG acoustic Courant condition on element size $\Delta x$. For high-density plasmas (such as RAM-C II at 61 km where $\omega_{p,\max} \approx 6.79 \times 10^{11}\,\mathrm{rad/s}$, $\tilde{\omega}_{p,\max} \approx 2266\,\mathrm{m}^{-1}$), set `solver_options.time_step` explicitly (e.g. `1e-4` to `1e-5` light-metres).
2. **Resolving shock gradients and skin depths:**
   Ensure element sizes inside the dispersive region resolve both the shortest wavelength of interest ($\ge 3$ order-4 elements or $\ge 5$ order-2 elements per wavelength) and sharp profile transitions such as the bow shock thickness $d_{\mathrm{sh}} \approx 0.32\,\mathrm{mm}$.
3. **Reference 1D reflection case:**
   See [`testData/testInputs/1D_Drude_RAMC_61km/1D_Drude_RAMC_61km.json`](../testData/testInputs/1D_Drude_RAMC_61km/1D_Drude_RAMC_61km.json) and [`scripts/ramc_61km_s11_skrf.py`](../scripts/ramc_61km_s11_skrf.py) for an end-to-end 1D $S_{11}$ verification against a 2400-layer `scikit-rf` transfer-matrix reference across $1\,\mathrm{MHz}$–$250\,\mathrm{GHz}$.
