# Physics and formulation: ADE CFS-CPML

## References

### Primary (implementation source of truth)

**Gedney, S. D. and Zhao, B.**, *An Auxiliary Differential Equation Formulation for the Complex-Frequency Shifted PML*, IEEE Trans. Antennas Propag., 58(3):838–847, 2010.  
DOI: [10.1109/TAP.2009.2037765](https://doi.org/10.1109/TAP.2009.2037765)

Key property for dgtd: ADE CFS-PML yields **complete first-order differential equations** in time. The **same time-marching scheme** used in the interior (RK4) applies to PML regions and auxiliary variables — no separate convolution update algorithm.

### Secondary (stretch profiles and intuition)

**Taflove**, *Computational Electrodynamics* — CPML section: complex stretch, grading, parameter tuning. Use for parameter sanity checks, **not** for spatial discretization (no Yee grid assumptions).

### Related (same ADE idea, different PDE)

Zhang & Shen, Geophysics 2010 — unsplit ADE CFS-PML for velocity–stress elastodynamics; cites Gedney & Zhao. Useful for understanding the **ADE split of 1/s** (Eq. 13–18 in that paper), which mirrors the Maxwell derivation structure.

---

## Conceptual picture

### Frequency domain (why PML exists)

Open domains require absorption at truncated boundaries. PML replaces real coordinates with **complex stretched coordinates** so outgoing waves decay exponentially inside a layer while the interface with the physical domain remains reflectionless (in the continuous model with ideal profiles).

For direction **d** (x, y, or z), the stretch is a complex function **s_d** of position and frequency.

### CFS vs classical PML

**Classical** stretch (schematic):

\[
s_d = \kappa_d + \frac{\sigma_d}{j\omega}
\]

**Complex-frequency-shifted (CFS)** stretch (Gedney / Taflove family):

\[
s_d = \kappa_d + \frac{\sigma_d}{\alpha_d + j\omega}
\]

| Parameter | Role |
|-----------|------|
| **κ_d ≥ 1** | Scaling / anisotropic slowing along stretched direction; improves evanescent and grazing absorption when κ > 1 |
| **σ_d ≥ 0** | Absorption; grows with depth into PML |
| **α_d ≥ 0** | Frequency shift; stabilizes late-time behavior; **α = 0** is classical PML limit |

At the **vacuum–PML interface**: **κ = 1**, **σ = 0**, **α = 0** → medium matches vacuum (reflectionless matching condition).

### Time domain (why ADE / CPML)

In frequency domain, stretch introduces rational functions of **jω**. Inverse Laplace transform gives **memory** — historically implemented as convolution integrals (**CPML** / recursive convolution).

**ADE formulation** introduces auxiliary memory variables **ψ** (paper uses **T** for stress components; Maxwell literature uses various ψ labels) with ODEs:

\[
\frac{\partial \psi}{\partial t} = -\alpha \psi + \beta \cdot (\text{field or spatial derivative coupling})
\]

The main field equations receive **correction terms** proportional to ψ and profile coefficients. This is **convolution PML without storing history**.

---

## ADE split (generic pattern from CFS literature)

Following the standard manipulation (Gedney & Zhao; Zhang & Shen Eq. 13–18 for analogous structure):

For a stretched derivative in direction **x**, write:

\[
\frac{1}{s_x} \frac{\partial}{\partial x} = \frac{1}{\kappa_x} \frac{\partial}{\partial x} + \text{(correction involving auxiliary } T_x \text{)}
\]

with auxiliary **T_x** satisfying in time domain (schematic):

\[
\frac{\partial T_x}{\partial t} = -\alpha_x T_x + \sigma_x \frac{\partial (\cdot)}{\partial x}
\]

(signs and which field component is differentiated must be taken from **Gedney & Zhao Maxwell ADE system**, not copied from elastodynamics).

The **correction** to the main PDE adds terms like **T_x / κ_x** (or equivalent) to the equation that would otherwise use the pure spatial derivative.

**Implementation rule:** When coding, open Gedney & Zhao and transcribe the **Maxwell curl** ADE system equation-for-equation. Do not guess signs from the seismic paper alone.

---

## Maxwell first-order system in dgtd (normalized vacuum)

Interior vacuum (no PML), normalized units (ε = μ = 1, c = 1):

\[
\frac{\partial \mathbf{E}}{\partial t} = \nabla \times \mathbf{H}, \qquad
\frac{\partial \mathbf{H}}{\partial t} = -\nabla \times \mathbf{E}
\]

(dgtd uses DG weak form; `GlobalEvolution` assembles curl + numerical flux via `DGOperatorFactory::buildGlobalOperator()`.)

**Conductive loss** (non-PML) today:

\[
\frac{\partial \mathbf{E}}{\partial t} \mathrel{+}= -\sigma_{\text{bulk}} \mathbf{E}
\]

via `collectGlobalConductiveOperator()` — **orthogonal** to PML stretch σ.

### PML-modified system (conceptual)

Inside PML-tagged elements, each active stretch direction **d** modifies how curl components involving **∂/∂x_d** enter the system:

1. **κ_d** part → can enter **inverse mass / material** side (effective ε, μ scaling along d).
2. **σ_d, α_d** part → **ADE auxiliaries** + correction terms in **`Mult()`** (PMLOperator).

Unsplit formulation: **same E, H components** as vacuum; no Berenger field splitting.

---

## Mapping Gedney ADE → dgtd operator split

| Continuous physics | dgtd implementation |
|--------------------|---------------------|
| Curl + DG flux (untouched by memory) | `globalOperator_` from `DGOperatorFactory` |
| κ stretch in constitutive / metric part | Modified **ε, μ** (or equivalent mass) on **PML element attributes** via `buildEpsMuPiecewiseVector()` |
| σ, α memory (ADE) | **`PMLOperator_`** or `applyPMLCoupling()` in `GlobalEvolution::Mult` |
| Auxiliary ψ evolution | Same **`Mult()`** output block for ψ DOFs |
| Vacuum–PML interface | Zero stretch at interface + standard interior DG flux |
| Outer PML boundary | **SMA** on domain boundary (Milestone A default) |

**Critical:** PML stretch **σ** must **not** be implemented as `Material::sigma_` / `buildSigmaPiecewiseVector()`.

---

## Auxiliary variable inventory (expectation)

Exact count follows Gedney & Zhao Maxwell formulation. Plan for:

- **Per active stretch direction** **d** at a point: auxiliaries coupling to derivatives **∂/∂x_d** appearing in curl.
- **Dimension-agnostic loop:** for each `d` in `{X,Y,Z}` with `d < mesh.Dimension()` and `d ∈ active_axes` for that material tag.
- **Corner elements** (e.g. tags with `"active_axes": ["X","Y"]`): combine profiles **κ_d(ρ_d), σ_d(ρ_d), α_d(ρ_d)** independently per direction; cross terms come from full Maxwell ADE, not ad-hoc corner code.

When **no PML tags** exist: **n_aux = 0**, no auxiliary blocks allocated.

---

## Profile functions (normalized coordinates)

Profiles are functions of **depth ρ_d** into the PML along direction **d**, measured from the **vacuum–PML interface** (see [03-json-schema-and-mesh.md](./03-json-schema-and-mesh.md)).

### Interface values (mandatory)

\[
\kappa_d(0) = 1, \quad \sigma_d(0) = 0, \quad \alpha_d(0) = 0
\]

### Depth normalization

For each PML element and direction **d**:

- **ρ_d = 0** at vacuum–PML interface face normal to **d**.
- **ρ_d = L_d** at the furthest point of that element along **+d** within the tagged region (or use quadrature-point geometry).

**L_d** comes from **mesh extent**, not JSON thickness.

### Grading (to be finalized at implementation from `target_reflection`)

JSON provides:

- `grading_order` (integer, e.g. 4)
- `kappa_max`, `alpha_max`
- `target_reflection` (e.g. 1e-6 → −60 dB design target; acceptance is −40 dB)

**Implementation task:** Derive **σ_max** (and κ profile if κ varies) from layer depth and `target_reflection` using Taflove/Gedney guidance, or use power law:

\[
\sigma_d(\rho) = \sigma_{d,\max} \left(\frac{\rho}{L_d}\right)^{m}, \quad m = \texttt{grading\_order}
\]

\[
\kappa_d(\rho) = 1 + (\kappa_{d,\max} - 1) \left(\frac{\rho}{L_d}\right)^{m}
\]

\[
\alpha_d(\rho) = \alpha_{d,\max} \left(\frac{\rho}{L_d}\right)^{m}
\]

**Verify against Gedney & Zhao** for exact profile conventions before locking constants.

---

## CPML vs UPML vs CFS-PML (terminology for Salvador)

| Term | Meaning in this project |
|------|-------------------------|
| **CFS-PML** | Stretch with **α + jω** in denominator |
| **CPML** | Time-domain via convolution **or equivalent ADE** |
| **ADE CFS-PML** | **What we implement** — Gedney & Zhao 2010 |
| **UPML** | Gedney uniaxial Maxwell form; same ADE system in unsplit variant |

Salvador suggested CPML + CFS-PML + Taflove: we implement **one** system: **ADE CFS-CPML** in DGTD.

---

## DGTD-specific notes (not FDTD)

1. **Spatial derivatives:** Weak DG derivatives + numerical fluxes, not finite differences on a staggered grid.
2. **High-order:** Polynomial order **p** on elements; PML profiles evaluated at **quadrature points** (dimension-agnostic).
3. **RK4 stages:** ψ must update inside **each** `Mult()` call at stage time — same as E, H (not SGBC checkpoint pattern).
4. **Curved elements:** Jacobian-aware derivatives already handled in `DerivativeIntegrator`; PML coefficients are **multipliers at QPs**, not a separate grid.

---

## Open items for implementer (check against paper)

- [x] Exact Maxwell ADE equation set (equation numbers in Gedney & Zhao 2010) — see **Gedney Maxwell ADE (dgtd transcription)** below
- [x] Sign convention for correction terms vs dgtd curl orientation — matched to `∂E/∂t = ∇×H`, `∂H/∂t = −∇×E`
- [x] Whether κ enters mass matrix only or also flux scaling in v1 — **Milestone A:** reference cases use `kappa_max=1`; κ-mass deferred; σ/α via ADE only
- [x] Formula linking `target_reflection` to σ_max for normalized units — implemented in `PMLProfiles.cpp`
- [x] Auxiliary DOF layout — **2 × ndofs × n_stretch_directions** (ψ_E and ψ_H per active axis); see `PMLAuxLayout`
- [x] Sign audit (2026-05-28): S1–S6 / global mult flip — **no t = 20 fix**; defaults unchanged — see [11-session-pml-sign-audit.md](./11-session-pml-sign-audit.md)

Document answers in PR or amend this file when resolved.

---

## Gedney Maxwell ADE (dgtd transcription)

Source: Gedney & Zhao, IEEE TAP 58(3):838–847, 2010 (DOI 10.1109/TAP.2009.2037765).  
Paper Eq. (1)–(3): CFS stretch \(s_d = \kappa_d + \sigma_d/(\alpha_d + j\omega)\).  
Paper Eq. (4)–(8): ADE split of \(1/s_d\) introducing auxiliary **b** (we label **ψ**).

### dgtd curl convention (matches `DGOperatorFactory`)

Normalized vacuum (\(\varepsilon=\mu=1\)):

\[
\frac{\partial \mathbf{E}}{\partial t} = \nabla \times \mathbf{H}, \qquad
\frac{\partial \mathbf{H}}{\partial t} = -\nabla \times \mathbf{E}
\]

Component form (\(d \in \{x,y,z\}\), skip \(d \ge \dim\)):

| Component | \((\nabla \times \mathbf{H})_c\) | \(-(\nabla \times \mathbf{E})_c\) |
|-----------|----------------------------------|-------------------------------------|
| \(E_x\) | \(\partial_y H_z - \partial_z H_y\) | — |
| \(E_y\) | \(\partial_z H_x - \partial_x H_z\) | — |
| \(E_z\) | \(\partial_x H_y - \partial_y H_x\) | — |
| \(H_x\) | — | \(\partial_y E_z - \partial_z E_y\) |
| \(H_y\) | — | \(\partial_z E_x - \partial_x E_z\) |
| \(H_z\) | — | \(\partial_x E_y - \partial_y E_x\) |

### ADE per active stretch direction \(d\)

For each direction \(d\) with active CFS stretch \((\kappa_d, \sigma_d, \alpha_d)\) at a point, introduce **two** auxiliary vectors collocated with DG DOFs:

- **ψ^E_d** — memory for terms where **∂E/∂x_d** appears in **Ḣ**
- **ψ^H_d** — memory for terms where **∂H/∂x_d** appears in **Ė**

**Auxiliary evolution** (Gedney ADE structure, paper Eq. (6)–(8) pattern):

\[
\frac{\partial \psi^E_d}{\partial t} = -\alpha_d\,\psi^E_d + \sigma_d\,\mathcal{D}_d(\mathbf{E})
\]

\[
\frac{\partial \psi^H_d}{\partial t} = -\alpha_d\,\psi^H_d + \sigma_d\,\mathcal{D}_d(\mathbf{H})
\]

where \(\mathcal{D}_d(\mathbf{F})\) is the **curl-coupled directional driver**: at each DOF, the linear combination of \(\partial F_c / \partial x_d\) with weights from the table above (only components whose curl row contains \(\partial/\partial x_d\)).

**Field corrections** (added to semidiscrete \(\dot{\mathbf{E}}\), \(\dot{\mathbf{H}}\) after `globalOperator_->Mult`):

\[
\frac{\partial E_c}{\partial t} \mathrel{+}= \sum_{d \in \text{active}} \frac{\psi^H_{d,c}}{\kappa_d}
\]

\[
\frac{\partial H_c}{\partial t} \mathrel{-}= \sum_{d \in \text{active}} \frac{\psi^E_{d,c}}{\kappa_d}
\]

(signs chosen so interface (\(\sigma=\alpha=0\), \(\kappa=1\)) gives zero correction).

**Spatial discretization:** \(\mathcal{D}_d\) uses the same weak directional derivative as `buildDerivativeSubOperator(d)` **plus** the matching interior one-normal face jump (see `collectGlobalOneNormalOperators`), both restricted to PML volume markers and scaled by \(\sigma\). Volume and face contributions enter with opposite sign on each field column (discrete SBP split). `PMLOperator_` is a preassembled CSR matrix built by `DGOperatorFactory::buildPMLOperator()` as `M_{\text{L2}}^{-1} \times` weak blocks (scalar inverse mass on ψ rows; Maxwell inverse mass on field correction rows), applied in `GlobalEvolution::applyPMLCoupling()` via `PMLOperator_->AddMult`.

### Auxiliary DOF count

\[
n_{\text{aux}} = 6 \times n_{\text{dofs}} \times n_{\text{stretch\_dirs}}
\]

Each stretch direction allocates **three** ψ^E blocks and **three** ψ^H blocks (one per vector component \(c \in \{x,y,z\}\)), each of length \(n_{\text{dofs}}\). Only component pairs with a nonzero curl entry (`pmlPsiEComponentActive` / `pmlPsiHComponentActive`) receive nonzero `PMLOperator_` blocks; uncoupled slots remain in state at zero.

where \(n_{\text{stretch\_dirs}}\) = number of distinct axes appearing in any region's `active_axes` (capped by `mesh.Dimension()`). Layout per stretch slot \(d\):

```text
[ Ex..Hz (6×N) | ψ^E_{d,X} (N) | ψ^E_{d,Y} (N) | ψ^E_{d,Z} (N) | ψ^H_{d,X} (N) | ψ^H_{d,Y} (N) | ψ^H_{d,Z} (N) | (next d …) ]
```

Example (1D TE, `active_axes: ["X"]`): only **ψ^E_{X,Z} ↔ Ey ↔ Hz** and **ψ^H_{X,Y} ↔ Hz ↔ Ey** are wired; **Hy** receives no Ey-driven ψ correction.

When no PML tags: \(n_{\text{aux}} = 0\).
